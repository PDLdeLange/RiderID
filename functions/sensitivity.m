function sen = sensitivity(sym,npm,mod)
% function [mod] = sensitivity(mod)
%
% input:  mod
% output: mod
%

    sel = logical(mod.X0);


    
    N = 2^12;
    D = 2^8;
    M = N/D;
    f = logspace(-2,2,N);
        
    X = mod.X;
    X1 = X(1);  X2 = X(2);  X3 = X(3);  X4 = X(4); %#ok<*NASGU>
    X5 = X(5);  X6 = X(6);  X7 = X(7);  X8 = X(8); X9 = X(9); 
    
    dKdXn = zeros([size(npm.y,2) sum(sel) N]);
    dydXn = zeros([size(npm.y,2) sum(sel) N]);
    tmpXn = zeros([sum(sel),sum(sel),size(npm.y,2),N]);    
    for k=1:N
       s = 1i*2*pi*f(k);
       symresults;
       dKdXn(:,:,k) = dKdXn_tmp(:,sel);
       dydXn(:,:,k) = dydXn_tmp(:,sel);
       disp([num2str(k) '/' num2str(N)]);
    end
    
    sen.f = f;
    sen.N = N;
    sen.M = M;
    sen.dKdXn = dKdXn;
    sen.dydXn = dydXn;
    
    % Multisine results
    N_bw = 8;
    f_bw = logspace(log10(0.2),log10(2),N_bw);
    Svv_bw = interp1(npm.f,npm.Svv,f_bw);

    % Determine covariance over specified frequency band
    dKdXn_bw = zeros([size(sym.dKdXn(:,sel)) N_bw]);
    dydXn_bw = zeros([size(sym.dydXn(:,sel)) N_bw]);
    tmp = NaN(sum(sel),sum(sel),2,N_bw);
    for k_bw = 1:N_bw;
        s = 1i*2*pi*f_bw(k_bw);
        symresults;
        dKdXn_bw(:,:,k_bw) = dKdXn_tmp(:,sel);
        dydXn_bw(:,:,k_bw) = dydXn_tmp(:,sel);
        for j = 1:2
            tmp(:,:,j,k_bw) =  1./Svv_bw(k_bw,j).*dydXn_bw(j,:,k_bw)'*dydXn_bw(j,:,k_bw); 
        end
    end
    cov = NaN(sum(sel),sum(sel),2);
    for j = 1:2
        cov(:,:,j) = 1/N_bw.*inv(sum(tmp(:,:,j,:),4)); % Parameter covariance
    end

    sen.cov_bw = cov;
    sen.N_bw = N_bw;
    sen.f_bw = f_bw;
    
%     sen.mf = mf;

end
