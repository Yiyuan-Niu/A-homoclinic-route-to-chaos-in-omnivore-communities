%% Figure 3e–f — Bifurcation Structure and Lyapunov Spectrum vs Control Parameter r
%
% This script explores the global dynamical behavior of the three-species
% model as the intrinsic growth rate r is varied continuously.
%
% For each value of r in r_list, the following steps are performed:
%   1) Integrate the system over a long time interval to remove transients.
%   2) Compute the full Lyapunov spectrum using the Benettin algorithm.
%   3) Extract local extrema (peaks) of R, C1, and C2 after transients.
%
% The results are used to construct:
%   • A Lyapunov exponent diagram (Fig. 3e), identifying chaotic regimes
%     via a positive maximal Lyapunov exponent.
%   • Classical bifurcation diagrams (Fig. 3f), obtained by plotting
%     asymptotic peaks of each state variable versus r.
%
%
    close all;
T_end = 2e4;

K0 = 10 ; 
w1 = 1 ; 
A1 = 0.23 ; 
B1 = 0.05 ; 
w2 = 0.05 ; 
A2 = 0.064 ; 
B2 = 0.25 ; 
w3 = 0.98 ; 
A3 = 0.14 ; 
B3 = 0.37 ; 
D1 = 0.3 ; 
D2 = 0.2 ; 

%% ------------------ Control parameter sweep ------------------
% Growth rate r is used as the bifurcation parameter
r_list = linspace(1.9, 3, 1000);

% Lyapunov spectrum
lambda_spectrum = zeros(length(r_list), 3);

%% ------------------ Storage for bifurcation data ------------------
% These arrays collect extrema (maxima/minima) of each variable
% after transient removal, to construct classical bifurcation diagrams
r_R = [];
R_all = [];

r_C1 = [];
C1_all = [];

r_C2 = [];
C2_all = [];

%% ================== Main parameter sweep loop ==================
% Parallel loop over r values for efficiency
parfor k = 1:length(r_list)
    %% Current parameter value
    r = r_list(k);

    args = [r, K0,...
    w1, A1, B1, ...
    w2, A2, B2, ...
    w3, A3, B3, ...
    D1, D2];

    %% Initial condition (kept identical for all r)
    begin = [6, 2, 13];

    sol_saddle = simulate_equation_1(args, begin,...
                        T_end, 'ode113'); 

    %% Lyapunov exponent computation (Benettin method)
    % Full Lyapunov spectrum is computed;
    [lambda_final, ~] = lyapunov_benettin(sol_saddle, args, 'ode113');

    lambda_spectrum(k, :) = lambda_final';

    disp(lambda_final)
    disp(k)

    %% Bifurcation
    % Only asymptotic dynamics are used for bifurcation analysis
    T_trans = 1e4;              
    idx = sol_saddle.t > T_trans;    
    t = sol_saddle.t(idx);
    R = sol_saddle.y(1, idx);
    C1 = sol_saddle.y(2, idx);  
    C2 = sol_saddle.y(3, idx);  

    % Extract local maxima and minima to build bifurcation diagrams
    R_extremum= findpeaks(R); 
    C1_extremum = findpeaks(C1); 
    C2_extremum = findpeaks(C2);

    if ~isempty(R_extremum)
        r_R = [r_R; r*ones(length(R_extremum),1)];
        R_all = [R_all; R_extremum(:)];
    end
    if ~isempty(C1_extremum)
        r_C1 = [r_C1; r*ones(length(C1_extremum),1)];
        C1_all = [C1_all; C1_extremum(:)];
    end
    if ~isempty(C2_extremum)
        r_C2 = [r_C2; r*ones(length(C2_extremum),1)];
        C2_all = [C2_all; C2_extremum(:)];
    end
    
end

%% ------------------ Lyapunov exponent diagram ------------------
% Fig. 3e
figure;
hold on
plot(r_list, lambda_spectrum(:,1), 'b-');
plot(r_list, lambda_spectrum(:,2), 'Color', [1 0.5 0]);
plot(r_list, lambda_spectrum(:,3), 'g-');
hold off
xlabel('r');
ylabel('lambda');
title('Bifurcation diagram: r - Lyapunov exponent');
grid on;


%% ------------------ Bifurcation diagram ------------------
figure;
scatter(r_R, R_all, 'filled','MarkerFaceColor',[1 0.5 0]);
xlabel('r');
ylabel('R');
title('Bifurcation diagram: r - Peak of R');
grid on;

% Fig. 3f: Bifurcation of C1
figure;
scatter(r_C1, C1_all, 'b');
xlabel('r');
ylabel('C1');
title('Bifurcation diagram: r - Peak of C1');
grid on;

% Fig. 3f: Bifurcation of C2
figure;
scatter(r_C2, C2_all, 'g');
xlabel('r');
ylabel('C2');
title('Bifurcation diagram: r - Peak of C2');
grid on;

%% Save Data
lambda_table = table(r_list(:), lambda_spectrum(:,1),lambda_spectrum(:,2),lambda_spectrum(:,3), ...
    'VariableNames', {'r', 'lambda_1', 'lambda_2', 'lambda_3'});
writetable(lambda_table, 'figure_3e.csv');

T_R = table(r_R(:), R_all(:), ...
    'VariableNames', {'r', 'R'});
writetable(T_R, 'figure_3f_R.csv');

T_C1 = table(r_C1(:), C1_all(:), ...
    'VariableNames', {'r', 'C1'});
writetable(T_C1, 'figure_3f_C1.csv');

T_C2 = table(r_C2(:), C2_all(:), ...
    'VariableNames', {'r', 'C2'});
writetable(T_C2, 'figure_3f_C2.csv');
