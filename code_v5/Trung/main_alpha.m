clear all; 
close all; 

%thêm cái này vào ko là banh xác, ko thêm thì ??ng h?i sao nó ko ch?y
addpath('C:\Users\ASUS\Dropbox\DATN-KyI-25_26\NHOMDATN\team_progress\Simulate\cvx');
cvx_setup
%thêm cái này vào, ko là banh xác, ko thêm thì ??ng h?i sao nó ko ch?y

snrdb = [20]; %dbm
ct=100; %to make the curve more smooth, this value should be increased
%M=2;
Mx=[2 :2: 10];
disc = 5; %BD users are within a square area [-disc disc].
d0=3;  %downlink user is located at [d0 0]
al=3;
alpha_vector =  [0.1 0.05 0.01 0.005 0.001];

sigma2_dbm= -94;%+10*log10(BW)+Nf; %Thermal noise in dBm -90+10+10
sigman=10^((sigma2_dbm-30)/10);


R0=3;%the target data rate is 2 bits/s/Hz   
eps0 = 2^R0-1;

for ki = 1:  length(alpha_vector)%length(Mx)
    M=5;%Mx(ki);
    alpha = alpha_vector(ki);
   Pmax = 10^((snrdb-30)/10);    
        
    sum1=0; sum2=0;  
  for ix = 1 :ct            
     %BD users
     hfading = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));
     locm = sign(rand(M,2)-0.5) .* disc.*rand(M,2); %location  of BD
     dist = max(1, sqrt(sum(locm.^2,2))); %distance to the base station
     hm = (abs(hfading)).^2./dist.^al; % channel for BD-BS
     
     gfading = complex(sqrt(0.5)*randn(M,1),sqrt(0.5)*randn(M,1));
     lxx = [d0*ones(M,1) zeros(M,1)];
     dist2 = max(1, sqrt(sum((locm-lxx).^2,2))); %distance between U0 and BD users
     gm = (abs(gfading)).^2./dist2.^al; % channel for BD-U0
     
     %downlink user
     h0fading = complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1));
     h0 = (abs(h0fading))^2/d0^al;
     
     %self interference 
     hsi = abs( complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1)) )^2;
     
     % signal
     s0= abs( complex(sqrt(0.5)*randn(1,1),sqrt(0.5)*randn(1,1)) )^2;
      
     %resource allocation the average one
     a = [0 ; hm.^2];
     
     d = [alpha*hsi; zeros(M,1)];
     a1 = eps0* (gm.*hm)';
     
     scale = d0^al;
     
     A = [-h0*scale a1*scale; zeros(M,1) -eye(M);-ones(M,1) eye(M); -1 zeros(1,M); 1 zeros(1,M)];
     b = [-eps0*sigman zeros(1,M) zeros(1,M) 0 Pmax ]';     

     %CVX part
     cvx_begin quiet
        variables y(M+1) z
        dual variables t1 t2 t3 
        minimize(  -a'*y )
        subject to
           t1 : A*y-b*z<=0;
           t3 : d'*y + sigman*z==1;
           t2 : z >= 0;
     cvx_end
    
     %rate one
     abar = [0 ; hm.^2*s0];
     p = y/z;
     ratex(ix) = log2(1+abar'*p/(d'*p +sigman));%log2(exp(1))*exp(1/(a'*y))*expint(1/(a'*y));
     
     
     %benchmark upper bound       
     %CVX part
     cvx_begin quiet
        variables y2(M+1) z2
        dual variables t21 t22 t23 
        minimize(  1 )
        subject to
           t21 : A*y2-b*z2<=0;
           t23 : d'*y2 + sigman*z2==1;
           t22 : z2 >= 0;
     cvx_end     
     p2 = y2/z2;
     sc=1;%rand;
     ratex2(ix) = log2(1+abar'*p2*sc/(d'*p2*sc +sigman));%log2(exp(1))*exp(1/(a'*y))*expint(1/(a'*y));
     
     
     %FD-OMA
     %resource allocation the average one
     ax = [0 ; hm(1).^2];     
     dx = [alpha*hsi; zeros(1,1)];
     ax1 = eps0* (gm(1).*hm(1));
     Ax = [-h0*scale ax1*scale; zeros(1,1) -eye(1);-ones(1,1) eye(1); -1 zeros(1,1); 1 zeros(1,1)];
     bx = [-eps0*sigman zeros(1,1) zeros(1,1) 0 Pmax ]';     

     %CVX part
     cvx_begin quiet
        variables y3(1+1) z3
        dual variables t31 t32 t33 
        minimize(  -ax'*y3 )
        subject to
           t31 : Ax*y3-bx*z3<=0;
           t33 : dx'*y3 + sigman*z3==1;
           t32 : z3 >= 0;
     cvx_end
     p3 = y3/z3;
     abarx = [0 ; hm(1)^2*s0];
     ratex3(ix) = log2(1+abarx'*p3/(dx'*p3 +sigman));%log2(exp(1))*exp(1/(a'*y))*expint(1/(a'*y));
  end 
  rate1(ki) = mean(ratex)
  rate2(ki) = mean(ratex2)
  rate3(ki) = mean(ratex3)
  
  
end


figure;

% S?A: Thay Mx b?ng alpha_vector
plot(alpha_vector, rate1, '-d', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(alpha_vector, rate2, '-s', 'LineWidth', 1.5, 'MarkerSize', 6);
plot(alpha_vector, rate3, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold off;

%% ===== LABEL ?ÚNG CHO FIG.5 =====
xlabel('FD Self-Interference Coefficient, \alpha', 'FontSize', 11);
ylabel('Average Sum Rate (bits/s/Hz)', 'FontSize', 11);
title('BAC-NOMA: Impact of FD Self-Interference (Fig.5)', 'FontSize', 12, 'FontWeight', 'bold');

legend('BAC-NOMA Proposed', 'BAC-NOMA Random', 'FD-OMA', 'Location', 'best');
grid on;

%% ===== THÊM SUBTITLE V?I THÔNG S? ?ÚNG =====
% S?A: ? trong subtitle ph?i là giá tr? hi?n t?i, ho?c ?? chung
subtitle_text = sprintf('M = %d BD devices, P_{max} = %d dBm, R_0 = %d BPCU', M, snrdb, R0);
subtitle(subtitle_text, 'FontSize', 10);

%% ===== ??O NG??C TR?C X (? càng nh? càng t?t) =====
set(gca, 'XDir', 'reverse');  % ? t? 0.1 ? 0.001 (l?n ??n nh?)
xlim([min(alpha_vector)*0.9, max(alpha_vector)*1.1]);

%% ===== ??NH D?NG TR?C X CHO ??P =====
% Chuy?n ? sang d?ng scientific notation ho?c gi? nguyên
xticks(alpha_vector);  % Hi?n th? ?úng các giá tr? ?
xticklabels(arrayfun(@(x) sprintf('%.3f', x), alpha_vector, 'UniformOutput', false));

% Ho?c dùng log scale n?u mu?n
% set(gca, 'XScale', 'log');

toc;

%% ===== IN K?T QU? RA COMMAND WINDOW =====
fprintf('\n========== K?T QU? FIG.5 (Rate vs ?) ==========\n');
fprintf('?\t\tFD-OMA\t\tRandom\t\tProposed\n');
fprintf('------------------------------------------------\n');
for i = 1:length(alpha_vector)
    fprintf('%.4f\t%.4f\t\t%.4f\t\t%.4f\n', ...
            alpha_vector(i), rate3(i), rate2(i), rate1(i));
end
fprintf('================================================\n');
 