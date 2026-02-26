function plot_actuators(ac_data)

    hold on; zoom on; grid on;

    % plot(ac_data.timestamp, ac_data.act_cmd_1);
    % plot(ac_data.timestamp, ac_data.act_state_1);
    % plot(ac_data.timestamp, ac_data.act_state_filt_1);
    
    plot(ac_data.timestamp, ac_data.act_state_1);
    plot(ac_data.timestamp, ac_data.act_state_2);
    plot(ac_data.timestamp, ac_data.act_state_3);
    plot(ac_data.timestamp, ac_data.act_state_4);
   
    xlabel('time [s]');
    ylabel('[pprz unit]');

    hold off;

end