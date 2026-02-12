function plot_time_series(sol, title_str)
    figure('WindowState','maximized');

    plot(sol.t, sol.y(1,:), 'Color',[1 0.5 0], 'DisplayName','R'); hold on;
    plot(sol.t, sol.y(2,:), 'b' ,'DisplayName','C1');
    plot(sol.t, sol.y(3,:), 'g', 'DisplayName','C2');

    xlim([1,max(sol.t)])
    set(gca,'XScale','log','YScale','log');
    xlabel('Time'); ylabel('Population');
    title(['Time series - ', title_str]);
    legend show;
end
