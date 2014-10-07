function [pR cRpd Rpd]=Rate_t_LMR(E, Jcoe, Reg, Redg, Reddg, Img, Imdg, Imddg, dtg, dtt)
    %INPUT:
    %E: exciton energy vector, unit in fs-1
    %Jcoe: coefficient matrix, dimension one represents the sites. 
    %Reg (lineshape real part, dt=0.1): each row contains the bath for each site.
    %Redg(lineshape first derivative real part)
    %Reddg(time correlation function realpart)
    %Img (lineshape imaginary part, dt=0.1): each row contains the bath for each site.
    %Imdg(lineshape first derivative imaginary part)
    %Imddg(time correlation function imaginary part) 
    
    %OUTPUT:
    %pR: a 3D structure of rate matrix predicted with modified Redfield theory, unit in fs-1
    %   dimension 1 and 2 are elements indicating transfer processes
    %   between exciton states; and dimension 3 are time, along which the
    %   rate are changed.
    %Rpd: the Markovian pure-dephasing rate matrix, dimension = 2.
    %cRpd:  the Non-Markovian pure-dephasing rate coefficient matrix, dimension = 2.
    %   The usage of cc is for time-dependent pure-dephasing rate Rpd(t).
    %   Each element of Rpd(t)=cRpd*Redg(t).
    %   All of the above matrice records EET rate from column sides to row
    %   sides. e.g. R(1,2) means from exciton 2 to exciton 1.

    %Rate integration:
    %unit of R sould be fs-1
    %======================================================================
    n = size(Jcoe,2);
    tmax=(size(Reg,2)-1)*dtg; %**maxima time calculated in bath
%     dt = tmax/(size(Redg,2)-1);
            
    L = -Imdg(:,end); %unit: fs-1, column vector
%     x = 0:dtg:tmax; %**time steps saved for each Reg, Img
    tt = 0:dtt:tmax; %**time steps saved for each Redg, Imdg
    
    nt=size(tt,2);

%     pReg=spline(x, Reg);
%     pRedg=spline(xx, Redg);
%     pReddg=spline(xx, Reddg);
% 
%     pImg=spline(x, Img);
%     pImdg=spline(xx, Imdg);
%     pImddg=spline(xx, Imddg);
    
    dG =Redg+1i*Imdg;
    
    %For cases that all sites use the same bath
    if size(Reg, 1)==1
        Jcoe=sum(Jcoe,1);
    else
    end
    
    dR = zeros(n,n,nt);
    cRpd = zeros(n,n);
    Rpd= zeros(n,n);
    
    for a = 1:n
        
        for b = [1:a-1 a+1:n]
            
            c = (2*Jcoe(:,b,b,a,a)-Jcoe(:,b,b,b,b)-Jcoe(:,a,a,a,a))'; %column to row
            va = (Jcoe(:,a,b,b,b)-Jcoe(:,a,b,a,a))';                  %column to row
            vb = (Jcoe(:,b,a,b,b)-Jcoe(:,b,a,a,a))';                  %column to row
            
            V = Jcoe(:,a,b,b,a)'*(Reddg+1i*Imddg)...
                    -(va.*vb)*dG.^2 ...
                    +1i*2*(((va.*Jcoe(:,b,a,a,a)'+vb.*Jcoe(:,a,b,a,a)')).*L')*dG...
                    +4*(Jcoe(:,b,a,a,a).*Jcoe(:,a,b,a,a))'*(L.^2);
                
     
            dR(b,a,:) = 2*real(exp(1i*(E(a)-E(b))*tt+c*Reg(1:dtt/dtg:end)+...
                1i*c*(L*tt+Img(1:dtt/dtg:end))).*V);
            cRpd(b,a) = -c;
            Rpd(b,a) = -c*Redg(:,end);         %Non-Markovian PD
            
% % %The Following code is the debugger:
% %             Rintgd=@(t) 2*real(exp(1i*(E(a)-E(b))*t).*AFv(t).*V(t));
% %             figure;
% %             plot(tt, Rintgd(tt),'b');
% %             hold on 
% %             plot(tt, imag(V(tt)), 'go',tt, real(V(tt)), 'g.');
% %             plot(tt, imag(AFv(tt)), 'ro',tt, real(AFv(tt)), 'r.');
% %             plot(tt, -imag(V(tt)).*imag(AFv(tt)), 'mo',tt, real(V(tt)).*real(AFv(tt)), 'm.');

        end
    end
    
    %Analytically integrate the rate equation==============================
    pR=spline(tt, dR);             %make a structure for pR
    clear dR
    
    pR.order = 5;                  %assign the order after integration
    coe = pR.coefs;                %save the coefficient for integration
    ndt = pR.pieces;               %ndt is the pieces;
    

    %integration
    pR.coefs = [coe(:,1)/4 coe(:,2)/3 coe(:,3)/2 coe(:,4) zeros(size(coe,1), 1)];
    coe = pR.coefs;
    coefs_sum = [coe(:,1)*dtt^4 coe(:,2)*dtt^3 coe(:,3)*dtt^2 coe(:,4)*dtt];
    pR.coefs = reshape(pR.coefs, n*n, ndt, 5);
    coefs_sum = reshape(coefs_sum, n*n, ndt, 4);
    pR.coefs(:, 2:ndt,5)= cumsum(sum(coefs_sum(:, 1:ndt-1, :), 3),2);
    pR.coefs = reshape(pR.coefs, n*n*ndt, 5);

end




