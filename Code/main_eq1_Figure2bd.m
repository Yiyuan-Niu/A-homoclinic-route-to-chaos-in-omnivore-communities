
%% =========================================================================
% Homoclinic Orbit in Three-Species Population Dynamics
% =========================================================================
% This script simulates a homoclinic orbit in a three-species ecological
% system, computes its Jacobian eigenvalues, and visualizes the 3D phase
% space trajectory for t > 1000 time units.
%    close all;
   
    r = 5.5;
    
    
    K0 = 10 ; 
    w1 = 0.7 ; 
    A1 = 0.35 ; 
    B1 = 0.2 ; 
    w2 = 0.1 ; 
    A2 = 0.45 ; 
    B2 = 0.35 ; 
    w3 = 0.53 ; 
    A3 = 0.67 ; 
    B3 = 0.64 ; 
    D1 = 0.35 ; 
    D2 = 0.55 ;
    
    

    args = [r, K0,...
    w1, A1, B1, ...
    w2, A2, B2, ...
    w3, A3, B3, ...
    D1, D2];

%% =========================================================================
% INITIAL CONDITIONS
% =========================================================================
    saddle_points = [ 0.7352,    2.5746,    1.4800];
    
    begin = zeros(1,3);
    begin(1) = saddle_points(1) * K0;
    begin(2) = saddle_points(2) / A1;
    begin(3) = saddle_points(3) / A2;

    begin


% Compute the Jacobian matrix at the saddle point and its eigenvalues.
    [eigvals, J] = compute_jac_eig_eq1(begin, args);
    eigvals

    sol_saddle = simulate_equation_1(args, begin,...
                        1e4, 'ode45');    

    plot_time_series(sol_saddle, 1);  
    
% Extract corresponding state variables, for t > 1e3
idx = sol_saddle.t > 1e3;

y = sol_saddle.y;
y = y(:, idx);
figure
plot3(y(1,:), y(2,:), y(3,:), 'LineWidth', 2)
grid on
xlabel('R')
ylabel('C_1')
zlabel('C_2')
title('Homoclinic orbit')

%% Save Data
sol_saddle_table = table(sol_saddle.t(:), ...
    sol_saddle.y(1,:)', ...
    sol_saddle.y(2,:)', ...
    sol_saddle.y(3,:)', ...
    'VariableNames', {'t', 'R', 'C1', 'C2'});

writetable(sol_saddle_table, 'figure_2b.csv');

Homoclinic = table(y(1,:)', ...
    y(2,:)', ...
    y(3,:)', ...
    'VariableNames', {'R', 'C1', 'C2'});

writetable(Homoclinic, 'figure_2d.csv');

