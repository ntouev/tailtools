function plot_ang_accel(ac_data)

    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.pdot_sp, LineWidth=1.5);
    h2 = plot(ac_data.timestamp, ac_data.pdot_filt, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('pdot [rad/s^2]');
    title('pdot');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h3 = plot(ac_data.timestamp, ac_data.qdot_sp, LineWidth=1.5);
    h4 = plot(ac_data.timestamp, ac_data.qdot_filt, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('qdot [rad/s^2]');
    title('qdot');


    ax3 = nexttile;
    hold on; zoom on; grid on;
    h5 = plot(ac_data.timestamp, ac_data.rdot_sp, LineWidth=1.5);
    h6 = plot(ac_data.timestamp, ac_data.rdot_filt, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('rdot [rad/s^2]');
    title('rdot');

    legend(ax1, [h1,h2], {'pdot sp','pdot'});
    legend(ax2, [h3,h4], {'qdot sp','qdot'});
    legend(ax3, [h5,h6], {'rdot sp','rdot'});

    linkaxes([ax1,ax2,ax3],'x');

end