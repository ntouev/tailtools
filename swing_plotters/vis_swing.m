function vis_swing(ac_data, fs_vis, anim)

    homeDir = getenv('HOME');
    addpath(genpath(fullfile(homeDir,'tailtools/stl')));

    vel_i = [ac_data.vel_n,ac_data.vel_e, ac_data.vel_d];
    norm_vel = sqrt(sum(vel_i.*vel_i, 2));

    ref = struct(); 
    actual = struct();

    ref.t = ac_data.timestamp;
    ref.p = [ac_data.pos_ref_n, ac_data.pos_ref_e, ac_data.pos_ref_d];
    ref.q = [ac_data.qs_sp, ac_data.qx_sp, ac_data.qy_sp, ac_data.qz_sp];
    
    actual.t = ac_data.timestamp;
    actual.p = [ac_data.pos_n, ac_data.pos_e, ac_data.pos_d];
    actual.q = [ac_data.qs, ac_data.qx, ac_data.qy, ac_data.qz];
    actual.AS = norm_vel;

    [ref, actual] = helper_interp(ref, actual, fs_vis);

    if anim == true
        fig1 = traj_vis_anim(ref, actual, 'traj', fs_vis, true);
    else
        fig1 = traj_vis_static(ref, actual, 'traj');
    end
    
    if anim == true
        set(fig1,'Position',[400, 100, 960, 1020]);
    else
        set(fig1,'Position',[60, 100, 800, 800]);
    end
end

%% helper
function [ref, actual] = helper_interp(ref, actual, fs_vis)
    t_vis = (ref.t(1) : 1/fs_vis : ref.t(end))';
    
    ref.p = [interp1(ref.t, ref.p(:,1), t_vis, 'linear','extrap'), ...
             interp1(ref.t, ref.p(:,2), t_vis, 'linear','extrap'), ...
             interp1(ref.t, ref.p(:,3), t_vis, 'linear','extrap')];
    ref.q = [interp1(ref.t, ref.q(:,1), t_vis, 'linear','extrap'), ...
             interp1(ref.t, ref.q(:,2), t_vis, 'linear','extrap'), ...
             interp1(ref.t, ref.q(:,3), t_vis, 'linear','extrap'), ...
             interp1(ref.t, ref.q(:,4), t_vis, 'linear','extrap')];
    ref.t = t_vis;
    
    actual.p = [interp1(actual.t, actual.p(:,1), t_vis, 'linear','extrap'), ...
                interp1(actual.t, actual.p(:,2), t_vis, 'linear','extrap'), ...
                interp1(actual.t, actual.p(:,3), t_vis, 'linear','extrap')];
    actual.q = [interp1(actual.t, actual.q(:,1), t_vis, 'linear','extrap'), ...
                interp1(actual.t, actual.q(:,2), t_vis, 'linear','extrap'), ...
                interp1(actual.t, actual.q(:,3), t_vis, 'linear','extrap'), ...
                interp1(actual.t, actual.q(:,4), t_vis, 'linear','extrap')];
    actual.AS = interp1(actual.t, actual.AS(:,1), t_vis, 'linear','extrap');
    actual.t = t_vis;
end