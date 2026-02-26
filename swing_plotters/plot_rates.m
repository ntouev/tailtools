function plot_rates(ac_data)

    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.rates_p_sp, LineWidth=1.5);
    h2 = plot(ac_data.timestamp, ac_data.rates_p, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('p [rad/s]');
    title('p');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h3 = plot(ac_data.timestamp, ac_data.rates_q_sp, LineWidth=1.5);
    h4 = plot(ac_data.timestamp, ac_data.rates_q, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('q [rad/s]');
    title('q');


    ax3 = nexttile;
    hold on; zoom on; grid on;
    h5 = plot(ac_data.timestamp, ac_data.rates_r_sp, LineWidth=1.5);
    h6 = plot(ac_data.timestamp, ac_data.rates_r, LineWidth=1.5);
    xlabel('time [s]');
    ylabel('r [rad/s]');
    title('r');

    legend(ax1, [h1,h2], {'p sp','p'});
    legend(ax2, [h3,h4], {'qsp','q'});
    legend(ax3, [h5,h6], {'r sp','r'});

    linkaxes([ax1,ax2,ax3],'x');

end