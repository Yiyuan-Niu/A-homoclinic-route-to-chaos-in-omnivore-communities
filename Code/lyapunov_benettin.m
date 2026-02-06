function [lambda_final, lambda_t] = lyapunov_benettin(sol, args, method)
% LYAPUNOV_BENETTIN  Compute Lyapunov exponents using the Benettin algorithm
%
% INPUTS:
%   sol    : structure containing the reference trajectory
%            sol.t : 1xN time vector
%            sol.y : 3xN state matrix [R; C1; C2]
%   args   : parameter for the system
%   method : ODE solver name, e.g., 'ode45', 'ode113', 'ode89'
% 
% OUTPUTS:
%   lambda_final : 3x1 vector of Lyapunov exponents at final time
%   lambda_t     : 3x(N-1) matrix of running Lyapunov exponents
% 
% ========================================================================
%               PERTURBED TRAJECTORY ( \tilde{y} ) using ode
%                  *--->----------------->*
%                 * |                      *
%                * /|                       *
%               * /                          *
% perturbation * /                            *
%   (delta_0) * /                              *    delta_t = \tilde{y} - y
%            * /                                *
%           * /                                  *
%          * /                                    *
%         *  eps0 * V_i                            *  eps0 * W_i
%        *                                          *
%       -->-----------------------------------------> t
%     t0          REFERENCE TRAJECTORY (y) using ode    t0+dt
% 
% ========================================================================

%   V_i  : 3 perturbation directions
%   W_i : Stetching of perturbations in 3 directions

%% Initialization
N = length(sol.t);          % number of time steps
lambda_t = zeros(3,N-1);    % running Lyapunov exponents
S = 0;                      % cumulative sum of logs

V = eye(3);                 % initial orthonormal perturbation directions

eps0 = 1e-8;                % initial perturbation magnitude

%% Loop over trajectory steps in sol
for idx = 1:N-1
    dt = sol.t(idx + 1) - sol.t(idx);

    % reference trajectory state at current step
    y0 = sol.y(:, idx);
    
    % store Stretching of perturbations in 3 Directions
    W = zeros(3,3);
    
    % Evolve perturbations along 3 directions
    for perturb_direction = 1:3

        % initial perturbation
        delta_0 = eps0 * V(:, perturb_direction);

        % perturbed initial condition \tilde{y0}
        yp = y0 + delta_0;
        
        % Compute perturbed trajector after dt
        sol_p = simulate_equation_1(args, yp', dt, method);
        delta_t = sol_p.y(:, end) - sol.y(:, idx + 1);
        
        % Normalize by initial perturbation magnitude
        % and store evolved Stretching of perturbations in Direction i
        W(:, perturb_direction) = delta_t / eps0;
    end
    
    % Re-orthonormalize perturbations using QR decomposition
    [Q, R] = qr(W); 
    V = Q;

    % Update cumulative sum of logarithms
    S = S + log(abs(diag(R)));
    
    % Compute Lyapunov exponents estimates under the current time node
    lambda_t(:,idx) = S / sol.t(idx+1);

end
      
    lambda_final = S/sol.t(end);
end
