clc;
clear all;

%% ================== SYSTEM PARAMETERS ==================
snrdb = 15:1:30;          % SNR in dB (15 -> 30 dB)
ct = 20;                % Monte-Carlo runs

M = 6;                    % Number of BD users
disc = 5;                 % User distribution range
d0 = 3;                   % Distance BS - DL user
al = 3;                   % Path-loss exponent

alpha = 0.01;             % Self-interference factor

sigma2_dbm = -94;         
sigman = 10^((sigma2_dbm-30)/10);   % Noise power (linear)

R0 = 3;                   
eps0 = 2^R0 - 1;          % SINR threshold

%% ================== RATE STORAGE ==================
rate_proposed = zeros(1,length(snrdb));
rate_upper    = zeros(1,length(snrdb));
rate_fdoma    = zeros(1,length(snrdb));

%% ================== MAIN LOOP ==================
for ki = 1:length(snrdb)

    Pmax = 10^((snrdb(ki)-30)/10);   % Max transmit power

    rate_mc  = zeros(1,ct);
    rate_mc2 = zeros(1,ct);
    rate_mc3 = zeros(1,ct);

    for ix = 1:ct

        %% ========== BACKSCATTER DEVICES ==========
        hfading = complex(randn(M,1),randn(M,1))/sqrt(2);
        locm = sign(rand(M,2)-0.5).*disc.*rand(M,2);
        dist = max(1, sqrt(sum(locm.^2,2)));
        hm = abs(hfading).^2 ./ dist.^al;

        gfading = complex(randn(M,1),randn(M,1))/sqrt(2);
        lxx = [d0*ones(M,1), zeros(M,1)];
        dist2 = max(1, sqrt(sum((locm-lxx).^2,2)));
        gm = abs(gfading).^2 ./ dist2.^al;

        %% ========== DOWNLINK USER ==========
        h0 = abs(complex(randn,randn)/sqrt(2))^2 / d0^al;

        %% ========== SELF INTERFERENCE ==========
        hsi = abs(complex(randn,randn)/sqrt(2))^2;

        %% ========== INFORMATION SIGNAL ==========
        s0 = abs(complex(randn,randn)/sqrt(2))^2;

        %% ================= PROPOSED SCHEME =================
        a  = [0; hm.^2];
        d  = [alpha*hsi; zeros(M,1)];
        a1 = eps0*(gm.*hm)';
        scale = d0^al;

        A = [-h0*scale      a1*scale;
              zeros(M,1)   -eye(M);
             -ones(M,1)     eye(M);
             -1             zeros(1,M);
              1             zeros(1,M)];

        b = [-eps0*sigman zeros(1,M) zeros(1,M) 0 Pmax]';

        cvx_begin quiet
            variables y(M+1) z
            maximize( a'*y )
            subject to
                A*y - b*z <= 0;
                d'*y + sigman*z == 1;
                z >= 0;
        cvx_end

        p = y/z;
        abar = [0; hm.^2*s0];
        rate_mc(ix) = log2(1 + (abar'*p)/(d'*p + sigman));

        %% ================= UPPER BOUND =================
        cvx_begin quiet
            variables y2(M+1) z2
            subject to
                A*y2 - b*z2 <= 0;
                d'*y2 + sigman*z2 == 1;
                z2 >= 0;
        cvx_end

        p2 = y2/z2;
        rate_mc2(ix) = log2(1 + (abar'*p2)/(d'*p2 + sigman));

        %% ================= FD-OMA (1 BD USER) =================
        ax  = [0; hm(1)^2];
        dx  = [alpha*hsi; 0];
        ax1 = eps0*gm(1)*hm(1);

        Ax = [-h0*scale  ax1*scale;
               0        -1;
              -1         1;
              -1         0;
               1         0];

        bx = [-eps0*sigman 0 0 0 Pmax]';

        cvx_begin quiet
            variables y3(2) z3
            maximize( ax'*y3 )
            subject to
                Ax*y3 - bx*z3 <= 0;
                dx'*y3 + sigman*z3 == 1;
                z3 >= 0;
        cvx_end

        p3 = y3/z3;
        abarx = [0; hm(1)^2*s0];
        rate_mc3(ix) = log2(1 + (abarx'*p3)/(dx'*p3 + sigman));

    end

    %% ========== MONTE-CARLO AVERAGE ==========
    rate_proposed(ki) = mean(rate_mc);
    rate_upper(ki)    = mean(rate_mc2);
    rate_fdoma(ki)    = mean(rate_mc3);
end

%% ================== LINEAR PLOT ==================
figure;
plot(snrdb, rate_fdoma,'-o','LineWidth',1.5); hold on;
plot(snrdb, rate_upper,'-s','LineWidth',1.5);
plot(snrdb, rate_proposed,'-^','LineWidth',1.5);
grid on;

xlabel('SNR (dB)');
ylabel('Average Sum Rate (bits/s/Hz)');
legend('FD-OMA','Upper Bound','Proposed','Location','NorthWest');
title('Average Sum Rate versus SNR (15–30 dB)');
