function plot_spec_force_cmd(ac_data)
    blue   = [0 0.4470 0.7410];
    orange = [0.8500 0.3250 0.0980];
    yellow = [0.9290 0.6940 0.1250];

    tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    norm_f = sqrt(ac_data.f_cmd_n.^2 + ac_data.f_cmd_e.^2 + ac_data.f_cmd_d.^2);

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.f_cmd_n, LineWidth=1.5);
    h2 = plot(ac_data.timestamp, ac_data.f_cmd_e, LineWidth=1.5);
    h3 = plot(ac_data.timestamp, ac_data.f_cmd_d, LineWidth=1.5);
    h4 = plot(ac_data.timestamp, norm_f, LineWidth=1);
    xlabel('time [s]');
    ylabel('f_i cmd [m/s^2]');
    title('Spec force cmd');

    legend(ax1, [h1,h2,h3,h4], {'fcmd_n','fcmd_e','fcmd_d','norm'});

    linkaxes([ax1],'x');

end