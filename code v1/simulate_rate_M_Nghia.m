function [Mx, rate_proposed, rate_upper, rate_fd_oma] = simulate_rate_vs_M(snrdb, ct, Mx, disc, d0, al, alpha_vec, sigma2_dbm, R0)
if nargin<1 || isempty(snrdb), snrdb = 20; end
if nargin<2 || isempty(ct), ct = 2000; end
if nargin<3 || isempty(Mx), Mx = 2:2:10; end
if nargin<4 || isempty(disc), disc = 5; end
if nargin<5 || isempty(d0), d0 = 3; end
if nargin<6 || isempty(al), al = 3; end
if nargin<7 || isempty(alpha_vec), alpha_vec = [0.1 0.05 0.01 0.005 0.001]; end
if nargin<8 || isempty(sigma2_dbm), sigma2_dbm = -94; end
if nargin<9 || isempty(R0), R0 = 3; end

% We'll pick a single alpha (matching original: alpha_vector(3) -> 0.01)
alpha = alpha_vec(min(3,length(alpha_vec)));

sigman = 10^((sigma2_dbm-30)/10);
eps0 = 2^R0 - 1;
Pmax = 10^((snrdb-30)/10);

% fraction of total Tx power reserved for downlink (heuristic)
p0_frac = 0.1;
pBD_frac = 1 - p0_frac;

% prealloc
rate_proposed = zeros(size(Mx));
rate_upper    = zeros(size(Mx));
rate_fd_oma   = zeros(size(Mx));

stream = RandStream('mt19937ar','Seed',sum(100*clock));


for k = 1:length(Mx)
    M = Mx(k);
    r_prop = zeros(ct,1);
    r_up   = zeros(ct,1);
    r_oma  = zeros(ct,1);
    for ix = 1:ct
        % --- fast channel generation
        hfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1));
        locm = sign(rand(M,2)-0.5).*disc.*rand(M,2);
        dist = max(1, sqrt(sum(locm.^2,2)));
        hm = (abs(hfading)).^2 ./ (dist.^al);   % BD-BS channel gains (power)
        
        gfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1));
        lxx = [d0*ones(M,1) zeros(M,1)];
        dist2 = max(1, sqrt(sum((locm - lxx).^2,2)));
        gm = (abs(gfading)).^2 ./ (dist2.^al);  % BD-U0 channel gains
        
        % downlink channel & interference
        h0fading = complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1));
        h0 = (abs(h0fading))^2 / d0^al;
        hsi = abs(complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1)))^2;
        s0 = abs(complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1)))^2;
        
        % heuristic power budget
        p0 = Pmax * p0_frac;
        P_bd_total = Pmax * pBD_frac;
        

        
        % ---------- Proposed: allocate BD power proportional to hm.^2 ----------
        w = hm.^2 + 1e-12;
        p_bd = P_bd_total * ( w / sum(w) );
        % numerator = sum(hm.^2 * s0 .* p_bd)
        num_prop = sum( (hm.^2 .* s0) .* p_bd );
        den_prop = alpha*hsi * p0 + sigman;
        r_prop(ix) = log2( 1 + num_prop / den_prop );
        
        % ---------- Upper bound (heuristic): all BD power to best contributor ----------
        % choose BD maximizing its contribution hm.^2 * s0
        [~, idx_best] = max(hm.^2 .* s0);
        p_bd_up = zeros(M,1); p_bd_up(idx_best) = P_bd_total;
        num_up = sum( (hm.^2 .* s0) .* p_bd_up );
        den_up = alpha*hsi * p0 + sigman;
        r_up(ix) = log2(1 + num_up / den_up);
        
       
        % ---------- FD-OMA (heuristic): choose BD maximizing gm.*hm ----------
        [~, idx_oma] = max(gm .* hm);
        p_bd_oma = zeros(M,1); p_bd_oma(idx_oma) = P_bd_total;
        num_oma = (hm(idx_oma)^2 * s0) * p_bd_oma(idx_oma);
        den_oma = alpha*hsi * p0 + sigman;
        r_oma(ix) = log2(1 + num_oma / den_oma);
    end
    
    rate_proposed(k) = mean(r_prop);
    rate_upper(k)    = mean(r_up);
    rate_fd_oma(k)   = mean(r_oma);
    
    fprintf('M=%2d | proposed=%.4f | upper=%.4f | fd-oma=%.4f\n', M, rate_proposed(k), rate_upper(k), rate_fd_oma(k));
end

% Plot (vs M)
figure;
plot(Mx, rate_fd_oma, '-o', Mx, rate_upper, '-s', Mx, rate_proposed, '-d','LineWidth',1.4);
xlabel('Number of BD users M');
ylabel('Average achievable rate (bits/s/Hz)');
legend('FD-OMA (heuristic)','Upper bound (heuristic)','Proposed (heuristic)','Location','best');
grid on;

end
