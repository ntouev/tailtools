function plot_f_model(ac_data, C_X, C_Z, C_T_v)
    blue   = [0 0.4470 0.7410];
    orange = [0.8500 0.3250 0.0980];
    yellow = [0.9290 0.6940 0.1250];

    quat = [ac_data.qs, ac_data.qx, ac_data.qy, ac_data.qz];
    vi = [ac_data.vel_n, ac_data.vel_e, ac_data.vel_d];
    vb = quatrotate(quat, vi);
    norm_v = sqrt(sum(vi.*vi, 2));

    fbx = C_X*norm_v.*vb(:,1);
    fbz = C_Z*norm_v.*vb(:,3) + C_T_v * ac_data.voltage.^2.*(ac_data.act_state_1.^2 + ac_data.act_state_2.^2 + ...
                                                             ac_data.act_state_3.^2 + ac_data.act_state_4.^2);
    fb = [fbx zeros(length(fbx),1) fbz];
    fi = quatrotate(quatinv(quat), fb);

    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax1 = nexttile;
    hold on; zoom on; grid on;
    h1 = plot(ac_data.timestamp, fb(:,1), LineWidth=1.5); hold on;
    h2 = plot(ac_data.timestamp, fb(:,2), LineWidth=1.5);
    h3 = plot(ac_data.timestamp, fb(:,3), LineWidth=1.5);
    xlabel('time [s]');
    ylabel('f_b [m/s^2]');
    title('f_b model');

    ax2 = nexttile;
    hold on; zoom on; grid on;
    h4 = plot(ac_data.timestamp, fi(:,1), LineWidth=1.5); hold on;
    h5 = plot(ac_data.timestamp, fi(:,2), LineWidth=1.5);
    h6 = plot(ac_data.timestamp, fi(:,3), LineWidth=1.5);
    xlabel('time [s]');
    ylabel('f_i [m/s^2]');
    title('f_i model');

    legend(ax1, [h1,h2,h3], {'fb_x','fb_y','fb_z'});
    legend(ax2, [h4,h5,h6], {'fi_x','fi_y','fi_z'});

    linkaxes([ax1,ax2],'x');

end