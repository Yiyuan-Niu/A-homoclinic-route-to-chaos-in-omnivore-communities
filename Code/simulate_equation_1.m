function [sol] = simulate_equation_1(args, y0, time_max, method)
% SIMULATE_EQUATION_1  Simulate a three-species population dynamics model
%
%   sol = simulate_equation_1(args, y0, time_max, method)
%
% INPUTS:
%   args      : Parameter vector [R0, K0, w1, A1, B1, w2, A2, B2, w3, A3, B3, D1, D2]
%               Parameters define species growth, interactions, and mortality
%   y0        : Initial conditions vector [R, C1, C2]
%   time_max  : Maximum simulation time
%   method    : ODE solver method, supports 'ode45' (default), 'ode113', 'ode89'
%
% OUTPUTS:
%   sol       : Structure containing simulation results
%               sol.t : Time vector
%               sol.y : Solution matrix (each row corresponds to R, C1, C2)
%
% MODEL DESCRIPTION (three-species interaction):
%
%         ┌───────────┐
%         │     R     │
%         └───────────┘
%          |         |
%   pred1  v         v  pred2
%        ┌─────┐   ┌─────┐
%        │ C1  │-> │ C2  │
%        └─────┘   └─────┘
%
% Equations:
%   dR/dt  = R0*R*(1 - R/K0) - A1*R*C1/(B1*R+1) - A2*R*C2/(B2*R+1)
%   dC1/dt = w1*A1*R*C1/(B1*R+1) - A3*C1*C2/(B3*C1+1) - D1*C1
%   dC2/dt = w2*A2*R*C2/(B2*R+1) + w3*A3*C1*C2/(B3*C1+1) - D2*C2
%
% Variables:
%   R  : Resource population
%   C1 : Consumer population 1
%   C2 : Consumer population 2
%
% ---------------------- Unpack parameters ----------------------
    [R0,K0,w1,A1,B1,w2,A2,B2,w3,A3,B3,D1,D2] = unpack_args(args);

    tspan = [0 time_max];
    opts = odeset('RelTol',1e-13,'AbsTol',1e-13);

% ---------------------- Solve ODE ----------------------
    switch method
        case 'ode89'
            [t,y] = ode89(@dPop_dt, tspan, y0,opts);
        case 'ode113'
            [t,y] = ode113(@dPop_dt, tspan, y0,opts);
        otherwise 
            [t,y] = ode45(@dPop_dt, tspan, y0,opts);
    end
% Transpose solution to have each row correspond to a variable
    sol.t = t;
    sol.y = y';

% ---------------------- ODE system ----------------------

    function dydt = dPop_dt(~,y)
        R  = y(1); % Resource
        C1 = y(2); % Consumer 1
        C2 = y(3); % Consumer 2
        
        % Resource growth and predation
        dR  = (R0*R*(1-R/K0) ...
              - A1*R*C1/(B1*R+1) ...
              - A2*R*C2/(B2*R+1));

        % Consumer 1 dynamics
        dC1 = (w1*A1*R*C1/(B1*R+1) ...
              - A3*C1*C2/(B3*C1+1) ...
              - D1*C1);

        % Consumer 2 dynamics
        dC2 = (w2*A2*R*C2/(B2*R+1) ...
              + w3*A3*C1*C2/(B3*C1+1) ...
              - D2*C2);
        
        dydt = [dR; dC1; dC2];
    end
end
