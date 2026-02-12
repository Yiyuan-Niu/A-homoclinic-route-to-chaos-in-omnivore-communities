%###########################################################
% MATLAB Script for Simulating Population Dynamics and Comparing
% with Experimental Data
%
% Author: [Your Name]
% Date: [YYYY-MM-DD]
% Purpose: 
%   Simulate a three-species system with different productivity
%   capacities, compute time-series and steady-state relative abundances,
%   and compare simulation results with experimental measurements.
%
% Workflow Overview:
%
%       ┌──────────────────────┐
%       │ Define parameters:   │
%       │ - Productivity levels│
%       │ - Species parameters │
%       │ - Interaction weights│
%       └─────────┬────────────┘
%                 │
%                 ▼
%       ┌──────────────────────┐
%       │ Loop over productivity│
%       │ levels (Low/Medium/High)│
%       └─────────┬────────────┘
%                 │
%                 ▼
%       ┌──────────────────────┐
%       │ Simulate ODE system: │
%       │ R, C1, C2 dynamics   │
%       │ using simulate_equation_1 │
%       └─────────┬────────────┘
%                 │
%                 ▼
%       ┌─────────────────────────────┐
%       │ Plot time-series (R, C1, C2)│
%       │ Compute steady-state mean   │
%       │ and error bounds (low/high) │
%       └─────────┬──────────────────┘
%                 │
%                 ▼
%       ┌─────────────────────────────┐
%       │ Compare simulation vs.      │
%       │ experimental relative abundances│
%       │ (scatter plots)             │
%       └─────────────────────────────┘
%
% Outputs:
%   1. Time series plots for each productivity level.
%   2. Steady-state mean and error bars of species abundances.
%   3. Scatter plots comparing model relative abundances with experimental data.
%
% Notes:
% - `simulate_equation_1(args, y0, t_end, solver)` is a user-defined function
%   that integrates the ODE system.
% - Steady-state statistics are computed after transient period T_trans = 5e4.
% - Error bars represent min/max deviation from the mean in the limit cycle.
% - Logarithmic scales are used for both axes in plots to better visualize dynamics.
% - Marker styles and colors distinguish productivity levels and species.
%
%###########################################################
    close all;
% Labels for three productivity regimes used in the study
    productivity = {'Low', 'Medium', 'High'};
    productivity_value = [0.06, 0.126, 0.164];
    marker_list= {'>', '*', 'o'};
% Preallocate matrices to store relative abundance comparison
model_store = zeros(3,3);   % rows = R,C1,C2 ; cols = Low,Med,High
exp_store   = zeros(3,3);
BC_similarity = zeros(1,3);

% Experimental data for three productivity levels
% [R, C1, C2]
    exp_data_list = {10.^[-3.8162,-4.6945,-5.4633], ...
                     10.^[-3.3695,-4.6773,-4.7589], ...
                     10.^[-3.1189,-4.6677,-4.3636]};

% Loop over the three productivity regimes
    for i = 1:3
        fprintf('=== %s capacity ===\n',  productivity{i});

% Experimental data relative abundance
        exp_data = exp_data_list{i};
        exp_rel = exp_data ./ sum(exp_data);
        
        r = productivity_value(i);
        K0 = 10;
    
        w1 = 0.7; A1 = 0.15; 
        w2 = 0.1; A2 = 0.1; 
        w3 = 0.2; A3 = 2; 
        
        D1 = 0.12; D2 = 0.1; 
      
        B1 = A1/7;
    
        B2 = A2/6;
    
        B3 = A3 / 3;
    
        args = [r, K0,...
        w1, A1, B1, ...
        w2, A2, B2, ...
        w3, A3, B3, ...
        D1, D2];
       
        
        y0 = [1.3, 0.35, 0.3];

        % Run long-term simulation to ensure asymptotic dynamics are reached
        [sol] = simulate_equation_1(args, y0, 1e5, 'ode45');
            
        t = sol.t;
        R = sol.y(1,:);
        C1 = sol.y(2,:);
        C2 = sol.y(3,:);
    
    [avg_abundance, lower_err_abundance, upper_err_abundance, ...
        model_rel] = compute_rel_abundance(sol);
    

    % Figure preview  
     % Two-panel layout: (1) dynamics, (2) model-vs-experiment comparison
        figure('WindowState','maximized');
        tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
        
        % --- Time series ---
        nexttile(1)
        
        hold on;
        plot(t, R, 'Color', [1 0.5 0], 'LineWidth', 1);
        plot(t, C1, 'b-', 'LineWidth', 1);
        plot(t, C2, 'g-', 'LineWidth', 1);
        hold off;
        
        % --- Extend axis for Errorbars---
        % --- Base position (log scale safe) ---
        xlim_current = xlim;
        x_base = xlim_current(2) * 1.1;
        
        
        x_R  = x_base * 1.00;
        x_C1 = x_base * 1.15;
        x_C2 = x_base * 1.30;
        
        hold on

    % Plot mean ± oscillation range for each variable
        R_mean = avg_abundance(1);
        C1_mean= avg_abundance(2);
        C2_mean = avg_abundance(3);

        R_err_low = lower_err_abundance(1);
        C1_err_low = lower_err_abundance(2);
        C2_err_low = lower_err_abundance(3);

        R_err_high = upper_err_abundance(1);
        C1_err_high = upper_err_abundance(2);
        C2_err_high= upper_err_abundance(3);

        h1 = errorbar(x_R, R_mean, ...
                      R_err_low, R_err_high, ...
                      marker_list{i}, ...
                      'LineWidth', 2, ...
                      'CapSize', 10, ...
                      'Color',[1 0.5 0],...
                      'MarkerFaceColor',[1 0.5 0]);
        
        h2 = errorbar(x_C1, C1_mean, ...
                      C1_err_low, C1_err_high, ...
                      marker_list{i}, ...
                      'LineWidth', 2, ...
                      'CapSize', 10, ...
                      'Color','b',...
                      'MarkerFaceColor','b');
        
        h3 = errorbar(x_C2, C2_mean, ...
                      C2_err_low, C2_err_high, ...
                      marker_list{i}, ...
                      'LineWidth', 2, ...
                      'CapSize', 10, ...
                      'Color','g',...
                      'MarkerFaceColor','g');
        
        hold off
        
        current = xlim;
        xlim([10, current(2)])
        
        xscale('log')
        yscale('log')
        xlabel('Time')
        ylabel('R')
        title(sprintf('=== %s capacity ===', productivity{i}))
        grid on
        
        nexttile(2)
    
    % ---- Bray–Curtis similarity ----
        num = sum(abs(model_rel(:) - exp_rel(:)));
        den = sum(model_rel(:) + exp_rel(:));
        
        BC_dissim = num / den;
        BC_similarity(i) = 1 - BC_dissim;
        
            
        hold on;
        scatter(exp_rel(1), model_rel(1), 120,...
                'MarkerFaceColor', [1 0.5 0],...
                'MarkerEdgeColor', [1 0.5 0],...
                 'Marker', marker_list{i})
        
        scatter(exp_rel(2), model_rel(2), 120,...
                    'MarkerFaceColor', 'b', ...
                    'MarkerEdgeColor', 'b', ...
                    'Marker', marker_list{i})
        
        scatter(exp_rel(3), model_rel(3), 120,...
                    'MarkerFaceColor','g', ...
                    'MarkerEdgeColor','g', ...
                    'Marker', marker_list{i})
        
        plot([0.01 1],[0.01 1],'k--','LineWidth',1.5)
        hold off;
        
        xlim([0.01 1])
        ylim([0.01 1])
        
        xscale('log')
        yscale('log')
        
        xlabel('Experimental relative abundance')
        ylabel('Model relative abundance (mean)')

        title(sprintf('BC similarity = %.3f', BC_similarity(i)))

        grid on
        axis square
    
    
    %% Save Data

    time_evolution_table = table(t, R', C1', C2');
    
    fig_suffix = char('a' + (i-1));

    filename = sprintf('figure_4%s.csv', fig_suffix);
    writetable(time_evolution_table, filename);

    model_store(:,i) = model_rel;
    exp_store(:,i)   = exp_rel(:);
    end

%% =========================
% Construct Table for Fig.4e
% ==========================

% Add an empty spacer column (for visual separation in Excel)
blank_col = repmat({''},3,1);

Fig4e_table = table( ...
    model_store(:,1), exp_store(:,1), blank_col, ...
    model_store(:,2), exp_store(:,2), blank_col, ...
    model_store(:,3), exp_store(:,3), ...
    'VariableNames', {'Low_model','Low_exp',' ', ...
                      'Medium_model','Medium_exp','  ', ...
                      'High_model','High_exp'} ...
    );

% Add row names (species)
Fig4e_table.Properties.RowNames = {'R','C1','C2'};

% Create a new row for Bray–Curtis similarity
BC_row = table( ...
     BC_similarity(1), BC_similarity(1), {''}, ...
     BC_similarity(2), BC_similarity(2), {''}, ...
     BC_similarity(3), BC_similarity(3),...
    'VariableNames', Fig4e_table.Properties.VariableNames);

BC_row.Properties.RowNames = {'BrayCurtis'};
Fig4e_table = [Fig4e_table; BC_row];

Fig4e_table
writetable(Fig4e_table, 'figure_4e.csv', 'WriteRowNames', true);


