function plot_ref_traj(ac_data)
    blue   = [0 0.4470 0.7410];
    orange = [0.8500 0.3250 0.0980];
    yellow = [0.9290 0.6940 0.1250];

    tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.pos_ref_n, LineWidth=1.5, Color=blue, LineStyle="-");
    h2 = plot(ac_data.timestamp, ac_data.pos_ref_e, LineWidth=1.5, Color=orange, LineStyle="-");
    h3 = plot(ac_data.timestamp, ac_data.pos_ref_d, LineWidth=1.5, Color=yellow, LineStyle="-");
    ylim([-6, 6]);
    xlabel('time [s]');
    ylabel('pos [m]');
    title('pos');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h4 = plot(ac_data.timestamp, ac_data.vel_ref_n, LineWidth=1.5, Color=blue, LineStyle="-");
    h5 = plot(ac_data.timestamp, ac_data.vel_ref_e, LineWidth=1.5, Color=orange, LineStyle="-");
    h6 = plot(ac_data.timestamp, ac_data.vel_ref_d, LineWidth=1.5, Color=yellow, LineStyle="-");
    xlabel('time [s]');
    ylabel('velocity [m/s]');
    title('velocity');

    ax3 = nexttile;
    hold on; zoom on; grid on;
    h7 = plot(ac_data.timestamp, ac_data.acc_ref_n, LineWidth=1.5, Color=blue, LineStyle="-");
    h8 = plot(ac_data.timestamp, ac_data.acc_ref_e, LineWidth=1.5, Color=orange, LineStyle="-");
    h9 = plot(ac_data.timestamp, ac_data.acc_ref_d, LineWidth=1.5, Color=yellow, LineStyle="-");
    xlabel('time [s]');
    ylabel('Acceleration [m/s^2]');
    title('Acceleration');

    ax4 = nexttile;
    hold on; zoom on; grid on;
    h10 = plot(ac_data.timestamp, ac_data.rates_ref_p, LineWidth=1.5, Color=blue, LineStyle="-");
    h11 = plot(ac_data.timestamp, ac_data.rates_ref_q, LineWidth=1.5, Color=orange, LineStyle="-");
    h12 = plot(ac_data.timestamp, ac_data.rates_ref_r, LineWidth=1.5, Color=yellow, LineStyle="-");
    xlabel('time [s]');
    ylabel('Angular rate [rad/s]');
    title('Angular rate');


    legend(ax1, [h1,h2,h3], {'p_N ref', 'p_E ref', 'p_D ref'});
    legend(ax2, [h4,h5,h6], {'v_N ref', 'v_E ref', 'v_D ref'});
    legend(ax3, [h7,h8,h9], {'a_N ref', 'a_E ref', 'a_D ref'});
    legend(ax4, [h10,h11,h12], {'p ref', 'q ref', 'r ref'});

    linkaxes([ax1,ax2,ax3,ax4],'x');

end