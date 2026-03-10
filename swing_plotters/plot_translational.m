function plot_translational(ac_data)

    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.pos_x, LineWidth=1.5);
    h2 = plot(ac_data.timestamp, ac_data.pos_y, LineWidth=1.5);
    h3 = plot(ac_data.timestamp, ac_data.pos_z, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('pos [m]');
    title('pos');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h4 = plot(ac_data.timestamp, ac_data.vel_x, LineWidth=1.5);
    h5 = plot(ac_data.timestamp, ac_data.vel_y, LineWidth=1.5);
    h6 = plot(ac_data.timestamp, ac_data.vel_z, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('velocity [m/s]');
    title('velocity');


    legend(ax1, [h1,h2,h3], {'N','E','D'});
    legend(ax2, [h4,h5,h6], {'N','E','D'});

    linkaxes([ax1,ax2],'x');

end