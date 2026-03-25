function plot_actuators(ac_data)
    blue   = [0 0.4470 0.7410];
    orange = [0.8500 0.3250 0.0980];
    yellow = [0.9290 0.6940 0.1250];

    hold on; zoom on; grid on;

    % plot(ac_data.timestamp, ac_data.act_cmd_1);
    % plot(ac_data.timestamp, ac_data.act_state_1);
    % plot(ac_data.timestamp, ac_data.act_state_filt_1);
    
    plot(ac_data.timestamp, ac_data.act_state_1);
    plot(ac_data.timestamp, ac_data.act_state_2);
    plot(ac_data.timestamp, ac_data.act_state_3);
    plot(ac_data.timestamp, ac_data.act_state_4);
    yline(9600, LineStyle="-.", LineWidth=2, Color="k");
    plot(ac_data.timestamp, ac_data.throttle, LineStyle=":", Color="r");

    ylim([0 10000]);

    ylabel('[pprz unit]');
    xlabel('time [s]');
    legend('Act state 1', 'Act state 2', 'Act state 3', 'Act state 4', 'saturation', 'throttle')

    hold off;

end