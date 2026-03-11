function plot_translational(ac_data)

    quat = [ac_data.qs, ac_data.qx, ac_data.qy, ac_data.qz];
    vi = [ac_data.vel_n, ac_data.vel_e, ac_data.vel_d];
    vb = quatrotate(quat, vi);

    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.pos_n, LineWidth=1.5);
    h2 = plot(ac_data.timestamp, ac_data.pos_e, LineWidth=1.5);
    h3 = plot(ac_data.timestamp, ac_data.pos_d, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('pos [m]');
    title('pos');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h4 = plot(ac_data.timestamp, ac_data.vel_n, LineWidth=1.5);
    h5 = plot(ac_data.timestamp, ac_data.vel_e, LineWidth=1.5);
    h6 = plot(ac_data.timestamp, ac_data.vel_d, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('velocity [m/s]');
    title('velocity');

    ax3 = nexttile;
    hold on; zoom on; grid on;
    h7 = plot(ac_data.timestamp, vb(:,1), LineWidth=1.5);
    h8 = plot(ac_data.timestamp, vb(:,2), LineWidth=1.5);
    h9 = plot(ac_data.timestamp, vb(:,3), LineWidth=1.5);
    xlabel('time [s]');
    ylabel('Body velocity [m/s]');
    title('Body velocity');


    legend(ax1, [h1,h2,h3], {'p_N','p_E','p_D'});
    legend(ax2, [h4,h5,h6], {'v_N','v_E','v_D'});
    legend(ax3, [h7,h8,h9], {'v_x','v_y','v_z'});

    linkaxes([ax1,ax2,ax3],'x');

end