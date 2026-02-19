function plot_stab_attitude(ac_data)

    tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.ang_accel_cmd_pdot, LineWidth=1.5);
    h2 = plot(ac_data.timestamp, ac_data.ang_accel_pdot, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('p dot [rad/s^2]');
    title('p dot');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h3 = plot(ac_data.timestamp, ac_data.rate_cmd_p, LineWidth=1.5);
    h4 = plot(ac_data.timestamp, ac_data.rate_filt_p, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('p [rad/s]');
    title('p');

    ax3 = nexttile;
    hold on; zoom on; grid on;
    h5 = plot(ac_data.timestamp, ac_data.att_cmd_phi, LineWidth=1.5);
    h6 = plot(ac_data.timestamp, ac_data.att_phi, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('phi [rad]');
    title('phi');

    ax4 = nexttile;
    hold on; zoom on; grid on;
    h7 = plot(ac_data.timestamp, ac_data.ang_accel_cmd_qdot, LineWidth=1.5);
    h8 = plot(ac_data.timestamp, ac_data.ang_accel_qdot, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('q dot [rad/s^2]');
    title('q dot');

    ax5 = nexttile;
    hold on; zoom on; grid on;
    h9 = plot(ac_data.timestamp, ac_data.rate_cmd_q, LineWidth=1.5);
    h10 = plot(ac_data.timestamp, ac_data.rate_filt_q, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('q [rad/s]');
    title('q');

    ax6 = nexttile;
    hold on; zoom on; grid on;
    h11 = plot(ac_data.timestamp, ac_data.att_cmd_theta, LineWidth=1.5);
    h12 = plot(ac_data.timestamp, ac_data.att_theta, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('theta [rad]');
    title('theta');

    ax7 = nexttile;
    hold on; zoom on; grid on;
    h13 = plot(ac_data.timestamp, ac_data.ang_accel_cmd_rdot, LineWidth=1.5);
    h14 = plot(ac_data.timestamp, ac_data.ang_accel_rdot, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('r dot [rad/s^2]');
    title('r dot');

    ax8 = nexttile;
    hold on; zoom on; grid on;
    h15 = plot(ac_data.timestamp, ac_data.rate_cmd_r, LineWidth=1.5);
    h16 = plot(ac_data.timestamp, ac_data.rate_filt_r, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('r [rad/s]');
    title('r');

    ax9 = nexttile;
    hold on; zoom on; grid on;
    h17 = plot(ac_data.timestamp, ac_data.att_cmd_psi, LineWidth=1.5);
    h18 = plot(ac_data.timestamp, ac_data.att_psi, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('psi [rad]');
    title('psi');

    % flight modes
    legend(ax1, [h1,h2], {'pdot ref','pdot'});
    legend(ax2, [h3,h4], {'p ref','p'});
    legend(ax3, [h5,h6], {'phi ref','phi'});
    legend(ax4, [h7,h8], {'qdot ref','qdot'});
    legend(ax5, [h9,h10], {'q ref','q'});
    legend(ax6, [h11,h12], {'theta ref','theta'});
    legend(ax7, [h13,h14], {'rdot ref','rdot'});
    legend(ax8, [h15,h16], {'r ref','r'});
    legend(ax9, [h17,h18], {'psi ref','psi'});

    linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9],'x');

end