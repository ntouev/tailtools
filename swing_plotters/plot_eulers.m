function plot_eulers(ac_data)

    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    quat = [ac_data.qs, ac_data.qx, ac_data.qy, ac_data.qz];
    eulZXY = quat2eul(quat,'ZXY');

    quat_sp = [ac_data.qs_sp, ac_data.qx_sp, ac_data.qy_sp, ac_data.qz_sp];
    eulZXY_sp = quat2eul(quat_sp,'ZXY');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, eulZXY_sp(:,2), LineWidth=1.5);
    h2 = plot(ac_data.timestamp, eulZXY(:,2), LineWidth=1.5);
    xlabel('time [s]');
    ylabel('phi [rad]');
    title('phi');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h3 = plot(ac_data.timestamp, eulZXY_sp(:,3), LineWidth=1.5);
    h4 = plot(ac_data.timestamp, eulZXY(:,3), LineWidth=1.5);
    xlabel('time [s]');
    ylabel('theta [rad]');
    title('theta');


    ax3 = nexttile;
    hold on; zoom on; grid on;
    h5 = plot(ac_data.timestamp, eulZXY_sp(:,1), LineWidth=1.5);
    h6 = plot(ac_data.timestamp, eulZXY(:,1), LineWidth=1.5);
    xlabel('time [s]');
    ylabel('psi [rad]');
    title('psi');

    legend(ax1, [h1,h2], {'phi sp','phi'});
    legend(ax2, [h3,h4], {'theta sp','theta'});
    legend(ax3, [h5,h6], {'psi sp','psi'});

    linkaxes([ax1,ax2,ax3],'x');

end