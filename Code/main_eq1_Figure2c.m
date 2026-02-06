
%% Figure 2c: Multiple trajectories under identical parameters
% This script simulates multiple trajectories of the three-species model
% under the same parameter set but different initial conditions.
%
% Purpose:
%   - Illustrate the Attracting Limit-Cycle
%   - Generate data corresponding to Figure 2c in the manuscript
%
% Model variables:
%   R  : Resource population
%   C1 : Consumer 1 population
%   C2 : Consumer 2 population

%% Clean up (open figure only once)    
close all;
   
%% Parameter definition

    r = 4.3;
    
    
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

%% Trajectory 1:
    begin = [10, 8, 2.2];

    trajectory_1 = simulate_equation_1(args, begin,...
                        1e4, 'ode89');    

    
%% Trajectory 2: 
    begin = [10, 12, 3];

    trajectory_2 = simulate_equation_1(args, begin,...
                        1e4, 'ode89');    

%% Plot Attracting Limit Cycle

    figure
    hold on
    
    plot3(trajectory_1.y(1,:), ...
        trajectory_1.y(2,:), ...
        trajectory_1.y(3,:), 'LineWidth', 2, 'Color', 'magenta')
    
    plot3(trajectory_2.y(1,:), ...
        trajectory_2.y(2,:), ...
        trajectory_2.y(3,:), 'LineWidth', 2, 'Color', [0.5 0 0.5])

    hold off
    grid on
    xlabel('R')
    ylabel('C_1')
    zlabel('C_2')
    title('Limit Cycle orbit')
    view(3)
    
% ---------------------- Camera settings for consistent visualization ----------------------
    ax = gca;

    ax.CameraPosition  = [67.8984  -6.6289   5.4422];
    ax.CameraTarget    = [ 6.1731   9.0000   2.7000];
    ax.CameraUpVector  = [ 0        0        1     ];
    ax.CameraViewAngle = 8.7139;

%% Save trajectory data (Figure 2c, trajectory 1 and 2)
trajectory_1_table = table(trajectory_1.y(1,:)', ...
    trajectory_1.y(2,:)', ...
    trajectory_1.y(3,:)', ...
'VariableNames', {'R', 'C1', 'C2'});

writetable(trajectory_1_table, 'figure_2c_trajectory_1.csv');

trajectory_2_table = table(trajectory_2.y(1,:)', ...
    trajectory_2.y(2,:)', ...
    trajectory_2.y(3,:)', ...
'VariableNames', {'R', 'C1', 'C2'});

writetable(trajectory_2_table, 'figure_2c_trajectory_2.csv');