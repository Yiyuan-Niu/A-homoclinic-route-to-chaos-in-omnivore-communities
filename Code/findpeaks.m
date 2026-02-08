function pks = findpeaks(y)
% FINDPEAKS  Simple peak and trough detector (pure MATLAB implementation)
%
% This function detects local extrema (both maxima and minima) in a
% one-dimensional signal without relying on MATLAB's Signal Processing
% Toolbox. It is designed for bifurcation analysis, where the goal is to
% extract asymptotic extrema of state variables rather than precise
% peak timing or prominence.
%
% INPUT:
%   y : 1D vector representing a time series (e.g., R(t), C1(t), or C2(t))
%
% OUTPUT:
%   pks : column vector containing unique local extrema values
%         (both local maxima and local minima)
%
% ALGORITHM:
%   • A point y(i) is classified as a local maximum if:
%         y(i) > y(i-1) AND y(i) > y(i+1)
%   • A point y(i) is classified as a local minimum if:
%         y(i) < y(i-1) AND y(i) < y(i+1)
%
%   All detected extrema are collected and then filtered using
%   uniquetol to remove nearly identical values caused by numerical
%   noise or finite precision.
%
% USAGE IN THIS PROJECT:
%   • Used to extract asymptotic extrema after transient removal
%   • Forms the basis of bifurcation diagrams (r vs peaks)
%
% NOTES:
%   • This implementation intentionally avoids thresholds, smoothing,
%     or prominence criteria to remain model-agnostic.
%   • Both maxima and minima are retained to capture full oscillation
%     envelopes in chaotic regimes.
%   • For long trajectories, this method is robust and computationally
%     inexpensive.
%

y = y(:); 
N = length(y);

pks = [];

%% Loop over trajectory steps in y
for i = 2:N-1
    if y(i) > y(i-1) && y(i) > y(i+1)
        pks = [pks; y(i)];
    end
    if y(i) < y(i-1) && y(i) < y(i+1)
        pks = [pks; y(i)];
    end
end
pks = uniquetol(pks,1e-5);

end
