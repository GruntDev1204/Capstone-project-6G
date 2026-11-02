clear all; 
close all; 

addpath('C:\Users\ASUS\Dropbox\DATN-KyI-25_26\NHOMDATN\team_progress\Simulate\cvx');
cvx_setup

tic;
snrdb = [20]; %dbm  (Gi? nguyên Pmax)
ct = 100; % S? l?n Monte Carlo

%base code
%Mx=[2 :4: 10];
%base code

M = 10;  % fixed M = 6
Mx = [10];  % Ch? 1 giá tr? M = 6

disc = 5; % BD users trong square area [-disc disc]
d0 = 3;  % downlink user t?i [d0 0]
al = 1;
alpha_vector = [0.1 0.05 0.01 0.005 0.001];
sigma2_dbm = -94;
sigman = 10^((sigma2_dbm-30)/10);
R0 = 3; % target data rate 3 bits/s/Hz   
eps0 = 2^R0 - 1;

for ki = 1:length(Mx)  % Gi? ch? ch?y 1 l?n (vì Mx ch? có [6])
    M = Mx(ki);  % M = 6
    alpha = alpha_vector(3);  % alpha = 0.01
    Pmax = 10^((snrdb-30)/10);    
        
    sum1 = 0; 
    sum2 = 0;  
    
    % Preallocate arrays
    ratex = zeros(1, ct);
    ratex2 = zeros(1, ct);
    ratex3 = zeros(1, ct);
    
    for ix = 1:ct            
        % BD users - T?O 6 DEVICES
        hfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1));
        locm = sign(rand(M,2)-0.5) .* disc .* rand(M,2); % location c?a 6 BD
        dist = max(1, sqrt(sum(locm.^2,2))); 
        hm = (abs(hfading)).^2 .* (dist.^(-al)); 
        
        % Channel BD-U0
        gfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1));
        lxx = [d0*ones(M,1) zeros(M,1)];
        dist2 = max(1, sqrt(sum((locm-lxx).^2,2))); 
        gm = (abs(gfading)).^2 .* (dist2.^(-al)); 
        
        % downlink user
        h0fading = complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1));
        h0 = (abs(h0fading))^2 * (d0^(-al));
        
        % self interference 
        hsi = abs(complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1)))^2;
        
        % signal
        s0 = abs(complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1)))^2;
       
        a = [0 ; hm.^2];  % 7x1 vector
        
        d = [alpha*hsi; zeros(M,1)];  % 7x1 vector
        a1 = eps0 * (gm.*hm)';  % 1x6 row vector
        
        scale = d0^al;
        
        % Xây d?ng A v?i ?úng kích th??c cho M=6
        A = [-h0*scale, a1*scale; 
             zeros(M,1), -eye(M);
             -ones(M,1), eye(M); 
             -1, zeros(1,M); 
             1, zeros(1,M)];
        
        b = [-eps0*sigman, zeros(1,M), zeros(1,M), 0, Pmax]';     
        
        % CVX part - BAC-NOMA Proposed
        cvx_begin quiet
            variables y(M+1) z  % y: 7x1, z: scalar
            dual variables t1 t2 t3 
            minimize(-a'*y)
            subject to
                t1: A*y - b*z <= 0;
                t3: d'*y + sigman*z == 1;
                t2: z >= 0;
        cvx_end
        
        % Tính rate
        abar = [0 ; hm.^2 * s0];  % ?? Ý: hm.^2 hay hm.^4?
        p = y/z;
        ratex(ix) = log2(1 + abar'*p / (d'*p + sigman));
        
        % Benchmark upper bound (Random)  
        cvx_begin quiet
            variables y2(M+1) z2
            dual variables t21 t22 t23 
            minimize(1)
            subject to
                t21: A*y2 - b*z2 <= 0;
                t23: d'*y2 + sigman*z2 == 1;
                t22: z2 >= 0;
        cvx_end     
        
        p2 = y2/z2;
        sc = 1; % rand n?u mu?n random
        ratex2(ix) = log2(1 + abar'*p2*sc / (d'*p2*sc + sigman));
        
        % FD-OMA (ch? dùng 1 BD device)
        % Ch?n device ??u tiên (ho?c random)
        m_oma = 1;  % Ch?n BD th? 1
        
        ax = [0 ; hm(m_oma).^2];     
        dx = [alpha*hsi; zeros(1,1)];
        ax1 = eps0 * (gm(m_oma).*hm(m_oma));
        Ax = [-h0*scale, ax1*scale; 
              zeros(1,1), -eye(1);
              -ones(1,1), eye(1); 
              -1, zeros(1,1); 
              1, zeros(1,1)];
        bx = [-eps0*sigman, zeros(1,1), zeros(1,1), 0, Pmax]';     

        cvx_begin quiet
            variables y3(2) z3  % 2 bi?n: P0 và P1
            dual variables t31 t32 t33 
            minimize(-ax'*y3)
            subject to
                t31: Ax*y3 - bx*z3 <= 0;
                t33: dx'*y3 + sigman*z3 == 1;
                t32: z3 >= 0;
        cvx_end
        
        p3 = y3/z3;
        abarx = [0 ; hm(m_oma)^2 * s0];
        ratex3(ix) = log2(1 + abarx'*p3 / (dx'*p3 + sigman));
    end 
    
    % Tính average rates
    rate1 = mean(ratex)   % BAC-NOMA Proposed
    rate2 = mean(ratex2)  % BAC-NOMA Random  
    rate3 = mean(ratex3)  % FD-OMA
    
    % ====== HI?N TH? K?T QU? CHO M=6 ======
    fprintf('\n========== K?T QU? M = %d ==========\n', M);
    fprintf('BAC-NOMA Proposed: %.4f bits/s/Hz\n', rate1);
    fprintf('BAC-NOMA Random:   %.4f bits/s/Hz\n', rate2);
    fprintf('FD-OMA:            %.4f bits/s/Hz\n', rate3);
    fprintf('==================================\n');
end

% ====== V? ?? TH? (ch? 1 ?i?m M=6) ======
figure;
bar([1, 2, 3], [rate3, rate2, rate1], 0.5);
set(gca, 'XTickLabel', {'FD-OMA', 'BAC-NOMA Random', 'BAC-NOMA Proposed'});
ylabel('Average Sum Rate (bits/s/Hz)');
title(sprintf('Performance M = %d BackCom Devices fixed', M));
grid on;

% Hi?n th? giá tr? trên các c?t
text(1, rate3, sprintf('%.3f', rate3), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(2, rate2, sprintf('%.3f', rate2), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(3, rate1, sprintf('%.3f', rate1), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

toc;