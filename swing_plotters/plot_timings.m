function plot_timings(ac_data)

    hold on; zoom on; grid on;

    plot(diff(ac_data.timestamp));
   
    xlabel('[cycle count]');
    ylabel('[sec]');

    hold off;

end