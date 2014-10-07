function [ t , Pn , RhoRec , Hp , disR , Gamma ] = EETinTetramer...
    ( r , theta , dt , Tmax , N )
% EETinTetramer calculates effective EET rate in tetramer
%   by means of quantum trajectory method
%   and Lindblad-form modified Redfield theory,
%   Note: Pure-dephasing process is indeed considered!
%
%   Definition of Input Parameters:
%   r: length between D1 and D2, in units of m
%   theta: (pi/2-theta) is angle between mu1 and center line
%   dt: time interval for evolution
%   Tmax: maximum evolution time
%   N: number of pieces in ensemble

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section for Inputing Parameters                             %
% in this section, initial condition is given                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

NH = 4; % dimension of Hilbert space

InitialState = [ 1; 0; 0; 0 ];
% initial state |psi(0)> in site basis [|1>;|2>;|3>;|4>]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section for Calculation                                     %
% in this section, calculation is performed                   %
% Caution: All energies and frequencies are in units of fs^-1 %
%          Time is in units of fs                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ H02 , pR , c_PD , pRedg ] = LMR_Para2(r,theta);
% H02: system Hamiltonian in site basis and in units of cm^-1
% pR: dissipation rate in matrix form
% c_PD & pRedg: used for calculation of pure-dephasing rate

[ vec , Hp ] = eig(H02);
% vec: exciton, eigen state in site basis
% Hs: system Hamiltonian in exciton basis
unitf = 2*pi*3e-5; % energy factor of cm^-1 to fs^-1
Hs = (Hp-1.2e4*eye(NH))*unitf; % eigen energy in units of fs^-1

State = vec'*InitialState; % wave function in exciton basis

Na = [ 0; 0; 0; 0; N ];
% population in different states: |1>,|2>,|3>,|4>, |psi(t)>
RhoRec = zeros(NH,NH,round(Tmax/dt));
% time evolution of density matrix in site basis

Pn = zeros(NH,round(Tmax/dt)); % time evolution of probabilities of excitons
%Nn = zeros(1,round(Tmax/dt)); % norm of determinstic wave function
R = zeros(NH); % Lindblad-form decay rate

t = dt:dt:Tmax;
Gamma = PureDephasingRate2(c_PD,pRedg,NH,t); % dephasing rate
disR = ppval(pR,t); % dissipation rate

% in this section, time evolution of density matrix of FMO is calculated
for t = dt : dt : Tmax
    
    tmpNa = Na; % temporary Na used for intermediate calculation
    count = t/dt;
    %Gamma = PureDephasingRate2(c_PD,pRedg,NH,t);
    %R = ppval(pR,t); % dissipation rate
    R = reshape(disR(:,:,count),NH,NH)+diag(Gamma(:,count)); 
    % decay rate including dissipation and dephasing
    totR = diag(sum(R,1)); % decay rate for nonHermitian Hamiltonian
    H = Hs-1i*totR/2; % total Hamiltonian for deterministic evolution
    
    % in this section, quantum jump for members in |psi(t)>
    if Na(NH+1,1) > 0
        rP = rand(Na(NH+1,1),1);
        Count = zeros(NH+1,1); % initialize target count
        
        % in this section,
        % quantum jump probability distribution is arranged
        dP = R.*(R>0).*repmat((abs(State).^2)',NH,1)*dt;
        dP = reshape(dP,NH^2,1);
        dP = 1-cumsum(dP);
        dP = [dP;0];
        
        % in this section, target state of quantum jump is determined
        [ii jj] = find(bsxfun(@plus,dP,-rP')<0);
        index = ii(logical([1;diff(jj)]));
        index = index(index<NH^2+1);
        target = rem(index,NH)+(rem(index,NH)==0)*NH;
        
        for j = 1 : NH
            Count(j,1) = CountNumber(target,j);
        end;
        tmpNa(1:NH,1) = tmpNa(1:NH,1)+Count(1:NH,1);
        tmpNa(NH+1,1) = tmpNa(NH+1,1)-sum(Count(:));
    end;
    if sum(tmpNa)~=N
        u=1;
    end;
    Rpc = R.*(R>0); % positive decay rate
    Rnc = -R.*(R<0); % negative decay rate
    
    % in this section, quantum jump for members in |k>
    for k = NH : -1 :1
        if Na(k,1) > 0
            rP = rand(Na(k,1),1);
            Count = zeros(NH+1,1); % initialize target count
            
            % in this section,
            % quantum jump probability distribution is arranged
            dPpc = Rpc(:,k)*dt;
            dPnc1NH = (Rnc(k,1:NH))'.*Na(1:NH,1)*dt/Na(k,1);
            dPnc0 = (Rnc(k,1:NH)').*(abs(State(1:NH,1)).^2)*dt*Na(NH+1,1)/Na(k,1);
            dP = 1-cumsum([ dPpc ; dPnc1NH ; dPnc0 ]);
            dP = [dP;0];
            [ii jj] = find(bsxfun(@plus,dP,-rP')<0);
            index = ii(logical([1;diff(jj)]));
            index = index(index<3*NH+1);
            
            target = (index<=NH).*index+(index>NH & index<=2*NH).*(NH-rem(index,NH))...
                +(index>2*NH)*(NH+1);
            for j = 1 : NH+1
                Count(j,1) = CountNumber(target,j);
            end;
            targetK = [1:k-1 k+1:NH+1];
            tmpNa(targetK,1) = tmpNa(targetK,1)+Count(targetK,1);
            tmpNa(k,1) = tmpNa(k,1)-sum(Count(targetK));
        end;
    end;
    if sum(tmpNa)~=N
        u=1;
    end;
    % wave function of deterministic evolution
    State = expm(-1i*H*dt)*State;
    State = State/norm(State);
    
    Na = tmpNa; % update state distribution due to quantum jump
    
    % in this section,
    % probabilities for different states are calculated
    count = round(t/dt);
    Pn(1 : NH,count) = (Na(1 : NH,1)+Na(NH+1,1)*(abs(State(1 : NH,1)).^2))/N;
    %Nn(count) = norm(State);
    RhoRec(:,:,count) = vec*(Na(NH+1,1)*(State*(State'))+diag(Na(1:NH,1)))*(vec')/N;
end;

t = dt:dt:Tmax;

toc;