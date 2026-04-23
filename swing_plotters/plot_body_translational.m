function plot_body_translational(ac_data)
    blue   = [0 0.4470 0.7410];
    orange = [0.8500 0.3250 0.0980];
    yellow = [0.9290 0.6940 0.1250];

    quat = [ac_data.qs, ac_data.qx, ac_data.qy, ac_data.qz];
    vi = [ac_data.vel_n, ac_data.vel_e, ac_data.vel_d];
    vb = quatrotate(quat, vi);
    norm_v = sqrt(sum(vi.*vi, 2));

    ai_filt = [ac_data.acc_filt_n, ac_data.acc_filt_e, ac_data.acc_filt_d];
    norm_a = sqrt(sum(ai_filt.*ai_filt, 2));

    fi_filt = ai_filt - [0 0 9.81];
    fb_filt = quatrotate(quat, fi_filt);
    norm_f = sqrt(sum(fb_filt.*fb_filt, 2));
    
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    
    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, vb(:,1), LineWidth=1.5);
    h2 = plot(ac_data.timestamp, vb(:,2), LineWidth=1.5);
    h3 = plot(ac_data.timestamp, vb(:,3), LineWidth=1.5);
    h4 = plot(ac_data.timestamp, norm_v, LineWidth=1, LineStyle='-');
    xlabel('time [s]');
    ylabel('Body velocity [m/s]');
    title('Body velocity');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h5 = plot(ac_data.timestamp, fb_filt(:,1), LineWidth=1.5);
    h6 = plot(ac_data.timestamp, fb_filt(:,2), LineWidth=1.5);
    h7 = plot(ac_data.timestamp, fb_filt(:,3), LineWidth=1.5);
    h8 = plot(ac_data.timestamp, norm_f, LineWidth=1, LineStyle='-');
    xlabel('time [s]');
    ylabel('f_b [m/s^2]');
    title('f_b filt');

    legend(ax1, [h1,h2,h3,h4], {'v_x','v_y','v_z','norm'});
    legend(ax2, [h5,h6,h7,h8], {'fb_x filt','fb_y filt','fb_z filt','norm'});

    linkaxes([ax1,ax2],'x');

end