function plot_translational(ac_data)
    blue   = [0 0.4470 0.7410];
    orange = [0.8500 0.3250 0.0980];
    yellow = [0.9290 0.6940 0.1250];

    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.pos_ref_n, LineWidth=1.5, Color=blue, LineStyle=":");
    h2 = plot(ac_data.timestamp, ac_data.pos_n, LineWidth=1.5, Color=blue);
    h3 = plot(ac_data.timestamp, ac_data.pos_ref_e, LineWidth=1.5, Color=orange, LineStyle=":");
    h4 = plot(ac_data.timestamp, ac_data.pos_e, LineWidth=1.5, Color=orange);
    h5 = plot(ac_data.timestamp, ac_data.pos_ref_d, LineWidth=1.5, Color=yellow, LineStyle=":");
    h6 = plot(ac_data.timestamp, ac_data.pos_d, LineWidth=1.5, Color=yellow);
    ylim([-6, 6]);
    xlabel('time [s]');
    ylabel('pos [m]');
    title('pos');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h7 = plot(ac_data.timestamp, ac_data.vel_sp_n, LineWidth=1.5, Color=blue, LineStyle=":");
    h8 = plot(ac_data.timestamp, ac_data.vel_n, LineWidth=1.5, Color=blue);
    h9 = plot(ac_data.timestamp, ac_data.vel_sp_e, LineWidth=1.5, Color=orange, LineStyle=":");
    h10 = plot(ac_data.timestamp, ac_data.vel_e, LineWidth=1.5, Color=orange);
    h11 = plot(ac_data.timestamp, ac_data.vel_sp_d, LineWidth=1.5, Color=yellow, LineStyle=":");
    h12 = plot(ac_data.timestamp, ac_data.vel_d, LineWidth=1.5, Color=yellow);
    xlabel('time [s]');
    ylabel('velocity [m/s]');
    title('velocity');

    % ax3 = nexttile;
    % hold on; zoom on; grid on;
    % h13 = plot(ac_data.timestamp, vb(:,1), LineWidth=1.5);
    % h14 = plot(ac_data.timestamp, vb(:,2), LineWidth=1.5);
    % h15 = plot(ac_data.timestamp, vb(:,3), LineWidth=1.5);
    % h16 = plot(ac_data.timestamp, norm_v, LineWidth=1.5);
    % xlabel('time [s]');
    % ylabel('Body velocity [m/s]');
    % title('Body velocity');

    ax4 = nexttile;
    hold on; zoom on; grid on;
    h17 = plot(ac_data.timestamp, ac_data.acc_sp_n, LineWidth=1.5, Color=blue, LineStyle=":");
    h18 = plot(ac_data.timestamp, ac_data.acc_filt_n, LineWidth=1.5, Color=blue);
    h19 = plot(ac_data.timestamp, ac_data.acc_sp_e, LineWidth=1.5, Color=orange, LineStyle=":");
    h20 = plot(ac_data.timestamp, ac_data.acc_filt_e, LineWidth=1.5, Color=orange);
    h21 = plot(ac_data.timestamp, ac_data.acc_sp_d, LineWidth=1.5, Color=yellow, LineStyle=":");
    h22 = plot(ac_data.timestamp, ac_data.acc_filt_d, LineWidth=1.5, Color=yellow);
    xlabel('time [s]');
    ylabel('Acceleration [m/s^2]');
    title('Acceleration');


    legend(ax1, [h1,h2,h3,h4,h5,h6], {'p_N ref', 'p_N', 'p_E ref', 'p_E', 'p_D ref', 'p_D'});
    legend(ax2, [h7,h8,h9,h10,h11,h12], {'v_N sp', 'v_N', 'v_E sp', 'v_E', 'v_D sp', 'v_D'});
    % legend(ax3, [h13,h14,h15,h16], {'v_x','v_y','v_z','norm'});
    legend(ax4, [h17,h18,h19,h20,h21,h22], {'a_N sp', 'a_N filt', 'a_E sp', 'a_E filt', 'a_D sp', 'a_D filt'});

    linkaxes([ax1,ax2,ax4],'x');

end