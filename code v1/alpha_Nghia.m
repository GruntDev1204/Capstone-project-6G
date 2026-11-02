
function [alpha_vector, rate_proposed, rate_upper, rate_fd_oma] = alpha(alpha_vector, ct, M, disc, d0, al, snrdb, sigma2_dbm, R0)
% ===== THAM S? M?C ??NH =====
% (Section IV - Numerical Studies trong bài báo)
if nargin < 1 || isempty(alpha_vector), alpha_vector = [0.1 0.05 0.01 0.005 0.001]; end  % Các giá tr? ? ?? v? Fig.5
if nargin < 2 || isempty(ct), ct = 100; end          % S? l??ng m?u Monte Carlo
if nargin < 3 || isempty(M), M = 6; end             % S? l??ng thi?t b? BackCom (BD) - Fig.4
if nargin < 4 || isempty(disc), disc = 5; end       % Vùng tri?n khai BD: [-5,5]m (Section IV)
if nargin < 5 || isempty(d0), d0 = 3; end           % Kho?ng cách BS t?i U0: 3m (Section IV)
if nargin < 6 || isempty(al), al = 3; end           % H? s? suy hao ???ng truy?n = 3 (Section IV)
if nargin < 7 || isempty(snrdb), snrdb = 30; end    % P_max = 20dBm (Section IV)
if nargin < 8 || isempty(sigma2_dbm), sigma2_dbm = -94; end  % Công su?t nhi?u = -94dBm (Section IV)
if nargin < 9 || isempty(R0), R0 = 3; end           % Rate yêu c?u t?i thi?u c?a U0 = 3 BPCU (Fig.5)

% ===== CHUY?N ??I ??N V? =====
% Eq.(3),(4): ?² là công su?t nhi?u
sigman = 10^((sigma2_dbm-30)/10);      % Chuy?n t? dBm sang Watt

% Eq.(P1d): P_max là công su?t phát t?i ?a c?a BS
Pmax = 10^((snrdb-30)/10);  % Chuy?n t? dBm sang Watt

% ===== PHÂN CHIA CÔNG SU?T HEURISTIC =====
% ?ây là heuristic t? ??t, KHÔNG có trong bài báo!
% Bài báo t?i ?u c? P0 và ?_m, code này ch? t?i ?u công su?t cho BD
p0_frac = 0.1;              % 10% công su?t cho downlink (P0)
pBD_frac = 1 - p0_frac;     % 90% công su?t cho các BD (P_bd)

% ===== KH?I T?O BI?N L?U K?T QU? =====
rate_proposed = zeros(size(alpha_vector));  % Ph??ng pháp ?? xu?t (CVX)
rate_upper    = zeros(size(alpha_vector));  % Upper bound heuristic
rate_fd_oma   = zeros(size(alpha_vector));  % FD-OMA heuristic

% ===== VÒNG L?P CHÍNH THEO ? =====
for k = 1:length(alpha_vector)
    alpha = alpha_vector(k);  % L?y giá tr? ? hi?n t?i (h? s? t? can nhi?u FD)
    
    % Kh?i t?o bi?n t?m cho Monte Carlo
    r_prop = zeros(ct,1);  % Rate c?a ph??ng pháp ?? xu?t
    r_up   = zeros(ct,1);  % Rate c?a upper bound  
    r_oma  = zeros(ct,1);  % Rate c?a FD-OMA
    
    % ===== VÒNG L?P MONTE CARLO =====
    for ix = 1:ct  % L?p ct l?n (m?c ??nh 100)
        
        %% === 1. T?O KÊNH TRUY?N NG?U NHIÊN ===
        % T?t c? theo mô hình Rayleigh fading
        
        % Kênh t? BD t?i BS: h_m (Eq.1)
        hfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1)); % Fading ph?c
        locm = sign(rand(M,2)-0.5).*disc.*rand(M,2); % V? trí ng?u nhiên c?a BD
        dist = max(1, sqrt(sum(locm.^2,2)));        % Kho?ng cách BD-BS
        hm = (abs(hfading)).^2 ./ (dist.^al);       % |h_m|² = |fading|² / d^?
        
        % Kênh t? BD t?i U0: g_m (Eq.2)
        gfading = complex(sqrt(0.5)*randn(M,1), sqrt(0.5)*randn(M,1));
        lxx = [d0*ones(M,1) zeros(M,1)];           % V? trí U0 c? ??nh: (d0, 0)
        dist2 = max(1, sqrt(sum((locm - lxx).^2,2))); % Kho?ng cách BD-U0
        gm = (abs(gfading)).^2 ./ (dist2.^al);     % |g_m|²
        
        % Kênh t? BS t?i U0: h_0 (Eq.4)
        h0fading = complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1));
        h0 = (abs(h0fading))^2 / d0^al;            % |h_0|²
        
        % Kênh t? can nhi?u (FD self-interference): h_SI (Eq.1)
        hsi = abs(complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1)))^2;
        
        % Tín hi?u downlink: s_0 (Eq.3) - gi? s? |s_0|² = 1
        s0 = abs(complex(sqrt(0.5)*randn(1,1), sqrt(0.5)*randn(1,1)))^2;
        
        %% === 2. PHÂN B? CÔNG SU?T C? ??NH ===
        % KHÁC v?i bài báo: Bài báo t?i ?u P0, code này c? ??nh P0 = 10% P_max
        p0 = Pmax * p0_frac;                % Công su?t downlink (P0)
        P_bd_total = Pmax * pBD_frac;      % T?ng công su?t cho các BD
        
        %% === 3. PH??NG PHÁP ?? XU?T (CVX OPTIMIZATION) ===
        % M?c tiêu: T?i ?a hóa rate uplink R_sum (Eq.3)
        % R_sum = log(1 + (? a_i * p_i) / denom)
        
        % Tính h? s? a_i = |h_m|? * |s_0|² (t? s? Eq.3)
        a = (hm.^2 .* s0);  % hm.^2 = |h_m|? vì hm ?ã là |h_m|²
        
        % Tính m?u s?: ?*P0*|h_SI|² + ?² (m?u s? Eq.3)
        denom = alpha*hsi * p0 + sigman;
        
        % ?n ??nh s? h?c: tránh chia cho 0
        small_eps = 1e-12;
        a = max(a, small_eps);
        denom = max(denom, small_eps);
        
        % T?i ?u l?i b?ng CVX (t??ng t? Problem P6 trong bài báo)
        p_bd_cvx = zeros(M,1);
        use_cvx = true;
        
        try
            cvx_begin quiet
                variable p_bd(M,1) nonnegative  % Bi?n c?n t?i ?u: công su?t các BD
                maximize( log( 1 + (a' * p_bd) / denom ) )  % T?i ?a hóa rate
                subject to
                    sum(p_bd) <= P_bd_total;  % Ràng bu?c t?ng công su?t
            cvx_end
            
            % Ki?m tra CVX có gi?i thành công không
            if ~exist('cvx_status','var') || ~(strcmpi(cvx_status,'Solved') || strcmpi(cvx_status,'Inaccurate/Solved'))
                use_cvx = false;
            else
                p_bd_cvx = p_bd;
            end
        catch
            use_cvx = false;  % CVX b? l?i
        end
        
        % Fallback n?u CVX th?t b?i
        if ~use_cvx
            w = hm.^2 + 1e-12;
            p_bd_cvx = P_bd_total * ( w / sum(w) );  % Phân b? t? l? v?i |h_m|?
        end
        
        % Tính rate ??t ???c
        num_prop = sum( a .* p_bd_cvx );  % T? s?: ? a_i * p_i
        r_prop(ix) = log2( 1 + num_prop / denom );  % Rate uplink (bps/Hz)
        
        %% === 4. UPPER BOUND HEURISTIC ===
        % Phân b? toàn b? công su?t cho BD t?t nh?t (theo a_i)
        % ?ây là GI?I H?N TRÊN vì b? qua m?i ràng bu?c giao thoa
        [~, idx_best] = max(a);  % BD có a_i l?n nh?t
        p_bd_up = zeros(M,1);
        p_bd_up(idx_best) = P_bd_total;  % Cho h?t công su?t vào BD này
        
        num_up = sum( a .* p_bd_up );
        r_up(ix) = log2( 1 + num_up / denom );
        
        %% === 5. FD-OMA HEURISTIC ===
        % Ch?n 1 BD duy nh?t (OMA) theo tiêu chí hm² * gm
        % hm²: kênh lên BS, gm: giao thoa xu?ng U0
        [~, idx_best2] = max(hm.^2 .* gm);  % Cân b?ng gi?a rate và giao thoa
        p_bd_oma = zeros(M,1);
        p_bd_oma(idx_best2) = P_bd_total;  % OMA: ch? 1 BD ???c truy?n
        
        num_oma = (hm(idx_best2)^2 * s0) * p_bd_oma(idx_best2);
        r_oma(ix) = log2( 1 + num_oma / denom );
    end
    
    %% === 6. TÍNH TRUNG BÌNH MONTE CARLO ===
    rate_proposed(k) = mean(r_prop);
    rate_upper(k)    = mean(r_up);
    rate_fd_oma(k)   = mean(r_oma);
    
    % In k?t qu?
    fprintf('alpha=%.5g: proposed=%.4f  upper=%.4f  fd-oma=%.4f\n', ...
            alpha, rate_proposed(k), rate_upper(k), rate_fd_oma(k));
end

%% === 7. V? ?? TH? ===
figure;
plot(alpha_vector, rate_fd_oma, '-o', ...
     alpha_vector, rate_upper, '-s', ...
     alpha_vector, rate_proposed, '-d', 'LineWidth',1.3);
xlabel('\alpha (h? s? t? can nhi?u)');
ylabel('Rate trung bình (bps/Hz)');
legend('FD-OMA (heuristic)','Upper bound (heuristic)','Proposed (CVX)', 'Location','best');
grid on;

end