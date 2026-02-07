%% Figure 3: Lyapunov Exponents and Sensitivity to Perturbations
% This script reproduces all panels of Figure 3 in the manuscript.
% It demonstrates chaotic dynamics in the three-species model by:
%
% (a) Simulating a long reference trajectory in time,
% (b) Measuring the exponential divergence between two nearby trajectories,
% (c) Visualizing the resulting strange attractor in 3D phase space,
% (d) Computing Lyapunov exponents using the Benettin algorithm.
%

% The system consists of:
%   R  : resource
%   C1 : consumer 1
%   C2 : consumer 2
%
% A small perturbation is applied to the initial condition of C1,
% and the divergence between the perturbed and reference trajectories
% is monitored over time. Exponential growth of this separation,
% together with a positive maximal Lyapunov exponent, confirms
% deterministic chaos.
%
% All numerical data used to generate the figures are exported to CSV
% files to ensure full reproducibility of the results.
%
% Numerical methods:
%   - ODE integration: MATLAB ode45
%   - Lyapunov exponents: Benettin algorithm (trajectory-based)
%   - Optional validation: Jacobian-based linear propagator
%
%% Clean up and simulation parameters


close all;
T_end = 1e4;   % simulation end time

% ------------------ System parameters ------------------
r = 2.5 ; 
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
    

    args = [r, K0,...
    w1, A1, B1, ...
    w2, A2, B2, ...
    w3, A3, B3, ...
    D1, D2];
    % 
    % saddle_points = [0.6849;    0.6552;    0.8131];
    % 
    % begin = zeros(1,3);
    % begin(1) = saddle_points(1) * K0;
    % begin(2) = saddle_points(2) / A1;
    % begin(3) = saddle_points(3) / A2;

%% Initial conditions
% Reference trajectory
    begin = [6, 2, 13];
    sol_saddle = simulate_equation_1(args, begin,...
                        T_end, 'ode45');

    t1 = sol_saddle.t;      
    Y1 = sol_saddle.y;
    

% Perturbed trajectory (small perturbation in C1)
    begin_perturbation = [6, 2 + 10^-4, 13];
    sol_perturbation = simulate_equation_1(args, ...
                        begin_perturbation, T_end, 'ode45');

    t2 = sol_perturbation.t;      
    Y2 = sol_perturbation.y;    

%% Calculate Difference between : 
% the Reference (y1) and the Perturbed trajectory (y2)
% y2 - y1
% Interpolate to common time grid if needed
    t_common = intersect(t1, t2);  
    if numel(t1)==numel(t2) && all(t1==t2)
        t = t1;
        Y1i = Y1;
        Y2i = Y2;
    else
        t = t1;   
        Y1i = Y1;
        Y2i = interp1(t2, Y2.', t).'; 
    end
    
%% Compute difference between trajectories in Fig. 3b
    delta = abs((Y2i - Y1i));

%% Compute Lyapunov exponents using Benettin method
    [lambda_final, lambda_t] = lyapunov_benettin(sol_saddle, args, 'ode45');
    lambda_final

% Jacaobian matrix method: 
%  OPTION:     Uncomment lines below to check the numerical bias
    % [lambda_final_jacobian, lambda_t_jacobian] = lyapunov_benettin_jacobian_method(sol_saddle, args);
    % lambda_final_jacobian

    times = sol_saddle.t(2:end);


%% Plot
% Fig. 3a : Chaotic trajectory
figure;
plot(t1, Y1(1,:), 'Color', [1 0.5 0], 'LineWidth', 0.1,'DisplayName','R');
hold on;
plot(t1, Y1(2,:), 'b-', 'LineWidth', 0.1,'DisplayName','C1');
plot(t1, Y1(3,:), 'g-', 'LineWidth', 0.1, 'DisplayName','C2');
hold off;
xlabel('Time')
ylabel('Abundance')
grid on
set(gca, 'YScale', 'log')

% Fig. 3b : Delta_R, C1, C2
figure;
plot(t, delta(1,:), 'Color', [1 0.5 0], 'DisplayName','\delta R');
hold on;
plot(t, delta(2,:), 'b-', 'DisplayName','\delta C1');
plot(t, delta(3,:), 'g-', 'DisplayName','\delta C2');
ax = gca;
xl = xlim;
xlim([1 xl(2)]);
hold off;
xlabel('Time')
ylabel('Delta')
grid on
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')


% Fig. 3c : Strange attractor in 3D phase space
% Extract corresponding state variables, for t > 5e3
idx = t1 > 5e3;
y = Y1(:, idx);

x = y(1,:);  % color will follow R value to improve the visualization
y_axis = y(2,:); % C1
z = y(3,:);  % 2C

figure
surf([x; x], [y_axis; y_axis], [z; z], [x; x], ...  % color = x
     'EdgeColor', 'interp', ...   % interpolate colors along line
     'FaceColor', 'none', ...     % hide surface face
     'LineWidth', 2)

grid on
xlabel('R')
ylabel('C_1')
zlabel('C_2')
title('Strange Smale attractor')

% Fig. 3d : Lyapunov exponent
figure;
plot(times, lambda_t(1,:)', 'b-', 'LineWidth', 1.2); hold on;
plot(times, lambda_t(2,:)', 'Color', [1 0.5 0],'LineWidth', 1.0);
plot(times, lambda_t(3,:)', 'g-', 'LineWidth', 1.0);
yline(lambda_final(1),'--k',sprintf('final1=%.4g',lambda_final(1)));

ax = gca;
xl = xlim;
xlim([1 xl(2)]);
xlabel('t'); ylabel('Lyapunov exponents (running)');
xscale('log')
legend('\lambda_1','\lambda_2','\lambda_3','Location','best');
title('Benettin method: Lyapunov convergence');
grid on;


%% Save trajectory data (Figure 2c, trajectory 1 and 2)
% Fig. 3a
fig3a_table= table(t1, ...
    Y1(1,:)', ...
    Y1(2,:)', ...
    Y1(3,:)', ...
'VariableNames', {'t', 'R', 'C1', 'C2'});

writetable(fig3a_table, 'figure_3a.csv');

% Fig. 3b
fig3b_table= table(t, ...
    delta(1,:)', ...
    delta(2,:)', ...
    delta(3,:)', ...
'VariableNames', {'t', 'delta_R', 'delta_C1', 'delta_C2'});

writetable(fig3b_table, 'figure_3b.csv');

% Fig. 3c
fig3c_table= table(y(1,:)', ...
    y(2,:)', ...
    y(3,:)', ...
'VariableNames', {'R', 'C1', 'C2'});

writetable(fig3c_table, 'figure_3c.csv');

% Fig. 3d
lypunov_table = table(times, ...
    lambda_t(1,:)', lambda_t(2,:)', lambda_t(3,:)', ...
   'VariableNames',{'time','lambda1','lambda2','lambda3'});
writetable(lypunov_table, 'figure_3d.csv')


lypunov_max_table = table(lambda_final(1,:), ...
    lambda_final(2,:), lambda_final(3,:), ...
   'VariableNames',{'lambda1','lambda2','lambda3'});
writetable(lypunov_table, 'figure_3d_last_value.csv')


