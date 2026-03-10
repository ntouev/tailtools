function plot_actuators(ac_data)

    hold on; zoom on; grid on;

    % plot(ac_data.timestamp, ac_data.act_cmd_1);
    % plot(ac_data.timestamp, ac_data.act_state_1);
    % plot(ac_data.timestamp, ac_data.act_state_filt_1);
    
    plot(ac_data.timestamp, ac_data.act_state_1);
    plot(ac_data.timestamp, ac_data.act_state_2);
    plot(ac_data.timestamp, ac_data.act_state_3);
    plot(ac_data.timestamp, ac_data.act_state_4);

    plot(ac_data.timestamp, ac_data.throttle, LineStyle=":", Color="r");

    yline(9600, LineStyle="--", Color="r", LineWidth=2);
   
    xlabel('time [s]');
    ylabel('[pprz unit]');

    legend('Act state 1', 'Act state 2', 'Act state 3', 'Act state 4', 'throttle', 'saturation')

    hold off;

end