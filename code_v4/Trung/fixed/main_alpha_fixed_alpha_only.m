clear all; close all; 

addpath('C:\Users\ASUS\Dropbox\DATN-KyI-25_26\NHOMDATN\team_progress\Simulate\cvx');
cvx_setup

% ====== THAM S? ======
snrdb = 20; 
ct = 100; 
disc = 5; 
d0 = 3;  
al = 3;
sigma2_dbm = -94;
sigman = 10^((sigma2_dbm-30)/10);
R0 = 3;   
eps0 = 2^R0 - 1;

% ====== FIXED M = 6 ======
M = 6;  % FIXED

% ====== 2 ? ======
alpha_values = [0.01, 0.001];

% ====== L?U K?T QU? ======
results = zeros(2, 3);  % 2 ? × 3 methods

% ====== SIMULATION ======
for alpha_idx = 1:2
    alpha_fixed = alpha_values(alpha_idx);
    Pmax = 10^((snrdb-30)/10);    
        
    ratex = zeros(1, ct);
    ratex2 = zeros(1, ct);
    ratex3 = zeros(1, ct);
    
    for ix = 1:ct            
        % BD users
        hfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1));
        locm = sign(rand(M,2)-0.5) .* disc .* rand(M,2);
        dist = max(1, sqrt(sum(locm.^2,2)));
        hm = (abs(hfading)).^2 ./ dist.^al;
        
        gfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1));
        lxx = [d0*ones(M,1) zeros(M,1)];
        dist2 = max(1, sqrt(sum((locm-lxx).^2,2)));
        gm = (abs(gfading)).^2 ./ dist2.^al;
        
        % downlink user
        h0fading = complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1));
        h0 = (abs(h0fading))^2 / d0^al;
        
        % self interference 
        hsi = abs(complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1)))^2;
        
        % signal
        s0 = abs(complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1)))^2;
        
        % resource allocation
        a = [0 ; hm.^2];
        d = [alpha_fixed*hsi; zeros(M,1)];
        a1 = eps0 * (gm.*hm)';
        scale = d0^al;
        
        A = [-h0*scale, a1*scale; 
             zeros(M,1), -eye(M);
             -ones(M,1), eye(M); 
             -1, zeros(1,M); 
             1, zeros(1,M)];
        b = [-eps0*sigman, zeros(1,M), zeros(1,M), 0, Pmax]';     
        
        % CVX - Proposed
        cvx_begin quiet
            variables y(M+1) z
            minimize(-a'*y)
            subject to
                A*y - b*z <= 0;
                d'*y + sigman*z == 1;
                z >= 0;
        cvx_end
        
        abar = [0 ; hm.^2 * s0];
        p = y/z;
        ratex(ix) = log2(1 + abar'*p / (d'*p + sigman));
        
        % Random
        cvx_begin quiet
            variables y2(M+1) z2
            minimize(1)
            subject to
                A*y2 - b*z2 <= 0;
                d'*y2 + sigman*z2 == 1;
                z2 >= 0;
        cvx_end     
        
        p2 = y2/z2;
        ratex2(ix) = log2(1 + abar'*p2 / (d'*p2 + sigman));
        
        % FD-OMA
        m_oma = 1;
        
        ax = [0 ; hm(m_oma).^2];     
        dx = [alpha_fixed*hsi; zeros(1,1)];
        ax1 = eps0 * (gm(m_oma).*hm(m_oma));
        Ax = [-h0*scale, ax1*scale; 
              zeros(1,1), -eye(1);
              -ones(1,1), eye(1); 
              -1, zeros(1,1); 
              1, zeros(1,1)];
        bx = [-eps0*sigman, zeros(1,1), zeros(1,1), 0, Pmax]';     

        cvx_begin quiet
            variables y3(2) z3
            minimize(-ax'*y3)
            subject to
                Ax*y3 - bx*z3 <= 0;
                dx'*y3 + sigman*z3 == 1;
                z3 >= 0;
        cvx_end
        
        p3 = y3/z3;
        abarx = [0 ; hm(m_oma)^2 * s0];
        ratex3(ix) = log2(1 + abarx'*p3 / (dx'*p3 + sigman));
    end 
    
    results(alpha_idx, 1) = mean(ratex);   % Proposed
    results(alpha_idx, 2) = mean(ratex2);  % Random
    results(alpha_idx, 3) = mean(ratex3);  % FD-OMA
end

% ====== V? FIGURE ??N GI?N ======

figure('Position', [100, 100, 1000, 400]);  % R?ng h?n ?? ch?a 2 subplots

% Subplot 1: ? = 0.01
subplot(1, 2, 1);  % 1 hàng, 2 c?t, v? trí 1
bar(1:3, results(1, :));
set(gca, 'XTick', 1:3, 'XTickLabel', {'Proposed', 'Random', 'FD-OMA'});
ylabel('Average Sum Rate (bits/s/Hz)');
title(sprintf('alpha = 0.01 \nM = %d, Pmax = %d dBm', M, snrdb));
grid on;

% Thêm giá tr? lên c?t
for i = 1:3
    text(i, results(1, i), sprintf('%.4f', results(1, i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom');
end

% Subplot 2: ? = 0.001
subplot(1, 2, 2);  % 1 hàng, 2 c?t, v? trí 2
bar(1:3, results(2, :));
set(gca, 'XTick', 1:3, 'XTickLabel', {'Proposed', 'Random', 'FD-OMA'});
ylabel('Average Sum Rate (bits/s/Hz)');
title(sprintf('alpha = 0.001\nM = %d, Pmax = %d dBm', M, snrdb));
grid on;

% Thêm giá tr? lên c?t
for i = 1:3
    text(i, results(2, i), sprintf('%.4f', results(2, i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom');
end

% ====== IN K?T QU? ======
fprintf('\n========== K?T QU? M = %d ==========\n', M);
fprintf('?\t\tProposed\tRandom\t\tFD-OMA\n');
fprintf('----------------------------------------\n');
fprintf('0.010\t%.4f\t\t%.4f\t\t%.4f\n', results(1,1), results(1,2), results(1,3));
fprintf('0.001\t%.4f\t\t%.4f\t\t%.4f\n', results(2,1), results(2,2), results(2,3));
fprintf('========================================\n');

toc;