function [eigvals, J] = compute_jac_eig_eq1(x, args)
% COMPUTE_JAC_EIG_EQ1  Compute the Jacobian matrix and its eigenvalues
%                       for the three-species model at a given state x.
%
%   [eigvals, J] = compute_jac_eig_eq1(x, args)
%
% INPUTS:
%   x    : 3x1 vector [R; C1; C2], current state of the system
%   args : Parameter vector, expected to contain at least 10 elements
%          [r, a1, a2, a3, a4, a5, a6, a7, a8, a9]
%
% OUTPUTS:
%   eigvals : 3x1 vector, eigenvalues of the Jacobian
%   J       : 3x3 Jacobian matrix evaluated at x
%

% ---------------------- Input validation ----------------------
if numel(x) ~= 3
    error('x must be 3x1 vector.');
end
if numel(args) < 10
    error('args must contain 10 elements [r,a1,...,a9].');
end

% ---------------------- Unpack variables ----------------------
R  = x(1); C1 = x(2); C2 = x(3);

[r, K0,...
w1, A1, B1, ...
w2, A2, B2, ...
w3, A3, B3, ...
D1, D2] = unpack_args(args);


% ---------------------- Precompute denominators ----------------------
denom_1R  = 1 + B1 * R;  % Used in C1-R interaction
denom_2R  = 1 + B2 * R;  % Used in C2-R interaction
denom_2C1 = 1 + B3 * C1; % Used in C1-C2 interaction


% ---------------------- Jacobian matrix elements ----------------------
% ∂(dR)/∂R
J11 = r*(1 - 2*R / K0) ...
    - A1 * C1/(denom_1R^2) ...
    - A2 * C2/(denom_2R^2);
% ∂dR/∂C1
J12 = - A1 * R / denom_1R;
% ∂dR/∂C2
J13 = - A2 * R / denom_2R;

% ∂dC1/∂R
J21 = w1 * A1 * C1 / (denom_1R^2);
% ∂dC1/∂C1
J22 = w1 * A1 * R / denom_1R ...
    - A3 * C2 / (denom_2C1^2) - D1;
% ∂dC1/∂C2
J23 = - A2 * C1 / denom_2C1;

% ∂dC2/∂R
J31 = w2* A2 * C2 / (denom_2R^2);
% ∂dC2/∂C1
J32 = w3 * A3 * C2 / (denom_2C1^2);
% ∂dC2/∂C2
J33 = w2 * A2 * R / denom_2R ...
    + w3 * A3 * C1 / denom_2C1 ...
    - D2;

% ---------------------- Assemble Jacobian ----------------------
% Jacobian J structure:
%   [ ∂R/∂R    ∂R/∂C1    ∂R/∂C2
%     ∂C1/∂R   ∂C1/∂C1   ∂C1/∂C2
%     ∂C2/∂R   ∂C2/∂C1   ∂C2/∂C2 ]

J = [ J11, J12, J13;
      J21, J22, J23;
      J31, J32, J33 ];

% ---------------------- Compute eigenvalues ----------------------
eigvals = eig(J);

end
