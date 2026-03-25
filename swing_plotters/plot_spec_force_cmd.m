function plot_spec_force_cmd(ac_data)
    blue   = [0 0.4470 0.7410];
    orange = [0.8500 0.3250 0.0980];
    yellow = [0.9290 0.6940 0.1250];

    tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.f_cmd_n, LineWidth=1.5);
    h2 = plot(ac_data.timestamp, ac_data.f_cmd_e, LineWidth=1.5);
    h3 = plot(ac_data.timestamp, ac_data.f_cmd_d, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('f_i cmd [m/s^2]');
    title('Spec force cmd');

    legend(ax1, [h1,h2,h3], {'fcmd_n','fcmd_e','fcmd_d'});

    linkaxes([ax1],'x');

end