%% Three-Species Dynamics Simulation
% This script simulates the three-species population dynamics model
% plots the time series, and saves the
% results to a CSV file.
%
% Model variables:
%   R  : Resource population
%   C1 : Consumer 1 population
%   C2 : Consumer 2 population
%
% System parameters:
%   r, K0 : Resource growth parameters
%   w1,A1,B1 : Consumer 1 interactions
%   w2,A2,B2 : Consumer 2 interactions
%   w3,A3,B3 : C1-C2 interaction
%   D1,D2 : Mortality rates

%% Clean up
    

close all;

%% Parameter setup
r = 4.3;         % intrinsic growth rate of resource
K0 = 10;         % resource carrying capacity

% Consumer 1
w1 = 0.7 ;     A1 = 0.35 ;      B1 = 0.2 ; 

% Consumer 2
w2 = 0.1 ;     A2 = 0.45 ;     B2 = 0.35 ; 

% C1-C2 interaction
w3 = 0.53 ;     A3 = 0.67 ;     B3 = 0.64 ; 

% Mortality
D1 = 0.35 ;     D2 = 0.55 ;
    

    args = [r, K0,...
    w1, A1, B1, ...
    w2, A2, B2, ...
    w3, A3, B3, ...
    D1, D2];

%% Initial condition
% Given saddle point in Eq S1 (\tilde{R}, \tilde{C1}, \tilde{C2})
    saddle_points = [ 0.6166,    2.7321,    1.3425];

% Rescale to actual initial conditions (R, C1, C2)
begin = zeros(1,3);
begin(1) = saddle_points(1) * K0;
begin(2) = saddle_points(2) / A1;
begin(3) = saddle_points(3) / A2;
    
begin
%% Simulate dynamics using ODE solver
sol_saddle = simulate_equation_1(args, begin,...
                        1e4,'ode45');    

%% Plot time series
plot_time_series(sol_saddle, 1);  


%% Save Data
sol_saddle_table = table(sol_saddle.t(:), ...
    sol_saddle.y(1,:)', ...
    sol_saddle.y(2,:)', ...
    sol_saddle.y(3,:)', ...
    'VariableNames', {'t', 'R', 'C1', 'C2'});

writetable(sol_saddle_table, 'figure_2a.csv');

    