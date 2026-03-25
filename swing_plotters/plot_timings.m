function plot_timings(ac_data)

    hold on; zoom on; grid on;

    % plot(ac_data.timestamp(2:end), diff(ac_data.timestamp));
    plot(ac_data.timestamp, ac_data.Ts);
    plot(ac_data.timestamp, ac_data.dt);
    hold off;
    
    xlabel('t [sec]');
    ylabel('[sec]');
    
    legend('Ts', 'dt');

end