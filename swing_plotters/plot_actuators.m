function plot_actuators(ac_data)

    hold on; zoom on; grid on;

    plot(ac_data.timestamp, ac_data.act_cmd_TL);
    plot(ac_data.timestamp, ac_data.act_cmd_TR);
    plot(ac_data.timestamp, ac_data.act_cmd_BR);
    plot(ac_data.timestamp, ac_data.act_cmd_BL);
    
    hold off;

end