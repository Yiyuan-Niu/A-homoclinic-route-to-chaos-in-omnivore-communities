function [lambda_final, lambda_t] = lyapunov_benettin_jacobian_method(sol, args)
% LYAPUNOV_BENETTIN_JACOBIAN_METHOD
% Compute Lyapunov exponents using the Benettin algorithm
% with linear propagation via the Jacobian matrix.
%
% This approach approximates the evolution of small perturbations
% using the Jacobian matrix at each time step instead of explicit
% integration of perturbed trajectories. It is useful for:
%   - Verifying zero Lyapunov exponent shifts
%   - Fast estimation of exponents for long trajectories
%   - Comparing with full trajectory propagation (simulate_equation_1)
%
%
% INPUTS:
%   sol  : structure containing the reference trajectory
%          sol.t : 1xN time vector
%          sol.y : 3xN state matrix [R; C1; C2]
%   args : parameter vector for the system, used in compute_jac_eig_eq1
%
% OUTPUTS:
%   lambda_final : 3x1 vector of final Lyapunov exponents
%   lambda_t     : 3x(N-1) matrix of running Lyapunov exponents
% 
% ========================================================================
%  \tilde{y0} : In jacobian matrix method,
%       ^       we do not directly intergate \tilde{y_t} from \tilde{y0};
%      /|\        
%       |                                    
%       |                                     
%       |          Jacobian Propagation
%       | delta_0 ---------------------> delta_t = expm(J * dt) * delta_0
%       |                                          
%      -->-----------------------------------------> t
%       t0        REFERENCE TRAJECTORY (y)         t0+dt
% ========================================================================
% 
% V_i  : 3 initial orthonormal perturbation directions
% W_i  : stretched perturbations after dt
% QR   : re-orthonormalization at each step 
% log(|Rii|) : accumulates growth rate -> Lyapunov exponent
% 
% This illustrates how small perturbations are evolved using the
% local Jacobian, then re-orthonormalized to compute Lyapunov exponents.
% ========================================================================

N = length(sol.t);
lambda_t = zeros(3,N-1);
S = zeros(3,1);

% initial orthonormal perturbation directions
% perturbation delta_0 in the direction of V_i
V = eye(3);

%% Loop over trajectory steps in sol
for idx = 1:N-1
    dt = sol.t(idx+1) - sol.t(idx);
    y0 = sol.y(:, idx);
    
    % Propagate perturbations using Jacobian
    [~, J] = compute_jac_eig_eq1(y0, args);  % 3x3 Jacobian
    
    % Linear propagation of perturbations
    W = expm(J * dt) * V;                    
        
    % QR decomposition to re-orthonormalize
    [Q, R] = qr(W,0);
    V = Q;
    
    % Update cumulative sum of logarithms
    S = S + log(abs(diag(R)));
    
    % Compute Lyapunov exponents estimates under the current time node
    lambda_t(:,idx) = S / sol.t(idx+1);
end

lambda_final = S / sol.t(end);

end
