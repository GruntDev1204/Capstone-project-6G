clear all; close all; 

addpath('C:\Users\ASUS\Dropbox\DATN-KyI-25_26\NHOMDATN\team_progress\Simulate\cvx');
cvx_setup

tic;
% ====== THAM S? ======
snrdb = 10:5:30;  % Pmax t? 10?30 dBm, b??c 5 dBm
ct = 300; 
% ====== FIXED M = 6 ======
M = 6;  
disc = 5; 
d0 = 3;  
al = 3;
alpha_vector = [0.1 0.05 0.01 0.005 0.001];
fixed_alpha = 0.01;
sigma2_dbm = -94;
sigman = 10^((sigma2_dbm-30)/10);
R0 = 3;   
eps0 = 2^R0 - 1;

% ====== PREALLOCATE ======
rate1 = zeros(1, length(snrdb));
rate2 = zeros(1, length(snrdb));
rate3 = zeros(1, length(snrdb));

% ====== LOOP QUA Pmax ======
for snr_idx = 1:length(snrdb)
    Pmax_dBm = snrdb(snr_idx);
    Pmax = 10^((Pmax_dBm-30)/10);  % dBm ? Watt
    alpha = fixed_alpha;  % ? = 0.01 c? ??nh
    
    fprintf('?ang x? lý Pmax = %d dBm, M = %d...\n', Pmax_dBm, M);
    
    % Preallocate cho Monte Carlo
    ratex_temp = zeros(1, ct);
    ratex2_temp = zeros(1, ct);
    ratex3_temp = zeros(1, ct);
    
    for ix = 1:ct            
        % BD users - LUÔN 6 DEVICES
        hfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1));
        locm = sign(rand(M,2)-0.5) .* disc .* rand(M,2);
        dist = max(1, sqrt(sum(locm.^2,2)));
        hm = (abs(hfading)).^2 .* (dist.^(-al));
        
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
        
        % resource allocation
        a = [0 ; hm.^2];
        d = [alpha*hsi; zeros(M,1)];
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
        
        % rate calculation
        abar = [0 ; hm.^2 * s0];
        p = y/z;
        ratex_temp(ix) = log2(1 + abar'*p / (d'*p + sigman));
        
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
        ratex2_temp(ix) = log2(1 + abar'*p2 / (d'*p2 + sigman));
        
        % FD-OMA
        m_oma = 1;
        
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
            variables y3(2) z3
            minimize(-ax'*y3)
            subject to
                Ax*y3 - bx*z3 <= 0;
                dx'*y3 + sigman*z3 == 1;
                z3 >= 0;
        cvx_end
        
        p3 = y3/z3;
        abarx = [0 ; hm(m_oma)^2 * s0];
        ratex3_temp(ix) = log2(1 + abarx'*p3 / (dx'*p3 + sigman));
    end 
    
    % L?u k?t qu?
    rate1(snr_idx) = mean(ratex_temp);
    rate2(snr_idx) = mean(ratex2_temp);
    rate3(snr_idx) = mean(ratex3_temp);
end

% ====== V? ?? TH? ======
figure('Position', [100, 100, 800, 500]);

plot(snrdb, rate1, '-d', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BAC-NOMA Proposed');
hold on;
plot(snrdb, rate2, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'BAC-NOMA Random');
plot(snrdb, rate3, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'FD-OMA');
hold off;

xlabel('Maximum Transmit Power, P_{max} (dBm)', 'FontSize', 12);
ylabel('Average Sum Rate (bits/s/Hz)', 'FontSize', 12);
title(sprintf('Rate vs P_{max} (M = %d BD Devices)', M), 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

% Thêm thông s?
subtitle_text = sprintf('? = %.3f, R_0 = %d BPCU, d_0 = %d m', alpha, R0, d0);
annotation('textbox', [0.15, 0.02, 0.7, 0.05], ...
           'String', subtitle_text, ...
           'FontSize', 10, ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center');

% ====== IN K?T QU? ======
fprintf('\n========== M = %d ==========\n', M);
fprintf('Pmax (dBm)\tProposed\tRandom\t\tFD-OMA\n');
fprintf('---------------------------------------------\n');
for snr_idx = 1:length(snrdb)
    fprintf('%d\t\t%.4f\t\t%.4f\t\t%.4f\n', ...
            snrdb(snr_idx), rate1(snr_idx), rate2(snr_idx), rate3(snr_idx));
end
fprintf('=============================================\n');

toc;