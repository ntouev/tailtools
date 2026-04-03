function plot_rates(ac_data)
    blue   = [0 0.4470 0.7410];
    orange = [0.8500 0.3250 0.0980];
    yellow = [0.9290 0.6940 0.1250];

    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, ac_data.rates_ref_p, LineWidth=1.5, LineStyle=":", Color=blue);
    h2 = plot(ac_data.timestamp, ac_data.rates_p_sp, LineWidth=1.5, Color=blue);
    h3 = plot(ac_data.timestamp, ac_data.rates_p, LineWidth=1.5, Color=orange);
    xlabel('time [s]');
    ylabel('p [rad/s]');
    title('p');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h4 = plot(ac_data.timestamp, ac_data.rates_ref_q, LineWidth=1, LineStyle=":", Color=blue);
    h5 = plot(ac_data.timestamp, ac_data.rates_q_sp, LineWidth=1.5, Color=blue);
    h6 = plot(ac_data.timestamp, ac_data.rates_q, LineWidth=1.5, Color=orange);
    xlabel('time [s]');
    ylabel('q [rad/s]');
    title('q');


    ax3 = nexttile;
    hold on; zoom on; grid on;
    h7 = plot(ac_data.timestamp, ac_data.rates_ref_r, LineWidth=1, LineStyle=":", Color=blue);
    h8 = plot(ac_data.timestamp, ac_data.rates_r_sp, LineWidth=1.5, Color=blue);
    h9 = plot(ac_data.timestamp, ac_data.rates_r, LineWidth=1.5, Color=orange);
    xlabel('time [s]');
    ylabel('r [rad/s]');
    title('r'); 

    legend(ax1, [h1,h2,h3], {'p ref', 'p sp', 'p'});
    legend(ax2, [h4,h5,h6], {'q ref', 'qsp', 'q'});
    legend(ax3, [h7,h8,h9], {'r ref', 'r sp', 'r'});

    linkaxes([ax1,ax2,ax3],'x');

end