function [avg_abundance, ...
    lower_err_abundance, upper_err_abundance, ...
    avg_relative] = compute_rel_abundance(sol)
%======================================================================
% compute_rel_abundance
%======================================================================
% Computes time-averaged absolute abundance and relative abundance of the
% state variables R, C1, C2 from an ODE solution structure, 
% The averaging is performed using time-weighted trapezoidal
% integration to ensure correctness for non-uniform solver time steps.

%----------------------------------------------------------------------
% INPUT
%----------------------------------------------------------------------
% sol : structure returned by MATLAB ODE solvers (e.g., ode45), containing:
%       sol.t  - time vector (1 × N or N × 1)
%       sol.y  - state matrix (n_state × N)
%                Each row is a variable (e.g., [R, C1, C2]),
%                each column corresponds to a time point.
%
%----------------------------------------------------------------------
% OUTPUT
%----------------------------------------------------------------------
% avg_abundance        : 1 × n_state vector
%                        Time-averaged absolute abundance of each variable
%                        over the post-transient window.
%
% lower_err_abundance  : 1 × n_state vector
%                        Deviation between the mean and the minimum value
%                        observed in the averaging window (useful for
%                        representing oscillation amplitude or uncertainty).
%
% upper_err_abundance  : 1 × n_state vector
%                        Deviation between the maximum value and the mean
%                        in the averaging window.
%
% avg_relative         : 1 × n_state vector
%                        Time-averaged relative abundance (composition)
%                        computed from normalized states:
%
%                            relative_i(t) = Yi(t) / sum_j Yj(t)
%
%----------------------------------------------------------------------
% METHOD
%----------------------------------------------------------------------
% 1. Convert solver output so each row corresponds to a time point:
%        Ys(t_k, :) = [R, C1, C2].
%
% 2. Compute instantaneous relative abundances by normalizing each row.
%
% 3. Select a post-transient averaging window:
%        Default: use all time points where t > 1e4.
%
%    This removes initial-condition effects and captures the asymptotic
%    regime (steady state or limit cycle).
%
% 4. If no samples satisfy t > 1e4 (e.g., short simulation), the function
%    automatically falls back to using the final 1% of simulated time.
%    This makes the routine robust to different simulation lengths.
%
% 5. Compute time-weighted averages using trapezoidal integration:
%
%        <Y> = ( ∫ Y(t) dt ) / (t_end − t_start)
%
%    This is critical because adaptive ODE solvers produce non-uniform
%    time spacing; a simple arithmetic mean would bias the result.
%
% 6. Estimate variability using extrema within the averaging window:
%
%        lower error = mean − min
%        upper error = max − mean
%
%    These values quantify oscillation amplitude for systems exhibiting
%    limit cycles rather than fixed equilibria.
%
%======================================================================
    Ys = sol.y';                  % each row = time point, columns = [R, C1, C2]
    
% Relative abundance
    total_abund = sum(Ys,2);
    rel = Ys ./ total_abund;
    
    tvec = sol.t(:);
    t_end = tvec(end);

% Select indices after a fixed transient cutoff (t > 10^4).
    idx_avg = find(tvec > 1e4);

    if isempty(idx_avg)
    % If the simulation did not extend beyond 1e4, fall back to a
    % relative-time criterion to remain robust to shorter simulations.

        warning('No time points with t > 1e4 found. Falling back to last 1%% of elapsed time.');
        t_start = tvec(1) + 0.99*(tvec(end) - tvec(1));  % begin of last 1% of elapsed time
        idx_avg = find(tvec >= t_start);
        % if still empty (extremely unlikely), use the last sample
        if isempty(idx_avg)
            warning('Last-1%% window also empty — using the last time sample.');
            idx_avg = length(tvec);
            t_start = tvec(idx_avg);
        else
            t_start = tvec(idx_avg(1));
        end
    else
        t_start = tvec(idx_avg(1));
    end


t_for_avg = tvec(idx_avg);
abundance_for_avg = Ys(idx_avg, :);   % m x 3
relative_for_avg = rel(idx_avg, :);


% -----------------------------
% Time-weighted averaging (critical for adaptive solvers)
% -----------------------------
% Preallocate vectors for averaged quantities.
avg_abundance = zeros(1, size(abundance_for_avg,2));
avg_relative = zeros(1, size(relative_for_avg,2));

for c = 1:size(abundance_for_avg,2)
    y = abundance_for_avg(:,c);
    avg_abundance(c) = trapz(t_for_avg, y) / (t_end - t_start);
    % Compute time-average using trapezoidal integration:
    %
    %   <Y_c> = ( ∫ Y_c(t) dt ) / Δt
    %
    % This avoids bias from non-uniform time steps generated by ode45.
end

lower_err_abundance = avg_abundance - min(abundance_for_avg);
upper_err_abundance = max(abundance_for_avg) - avg_abundance;

% Relative abundance trajectory of species R, C1, C2.
for c = 1:size(relative_for_avg,2)
    y = relative_for_avg(:,c);
    avg_relative(c) = trapz(t_for_avg, y) / (t_end - t_start);
end

end
