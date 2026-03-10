clear;

homeDir = getenv('HOME');
addpath(genpath(fullfile(homeDir,'tailtools/swing_plotters')));

meta = struct();

%%% real flight data
% ac_data = readtable('~/LOGS/swing/20260219_tuning/rpy.csv');
% ac_data = readtable('~/LOGS/swing/20260219_tuning/yaw_issue.csv');

% ac_data = readtable('~/LOGS/swing/20260304_stab/01-maiden.csv');
% ac_data = readtable('~/LOGS/swing/20260304_stab/02.Kq2_Komega10.csv');
% ac_data = readtable('~/LOGS/swing/20260304_stab/08.Kq2.5Komega12.csv');

% ac_data = readtable('~/LOGS/swing/20260306_more_stab/01.rigid_structure.csv');
% ac_data = readtable('~/LOGS/swing/20260306_more_stab/02.coupling_issue_1st_n_last_pitch.csv');
% ac_data = readtable('~/LOGS/swing/20260306_more_stab/03.training_dataset.csv');
% ac_data = readtable('~/LOGS/swing/20260306_more_stab/05.diff_eff_values.csv');

% ac_data = readtable('~/LOGS/swing/20260309_again_stab/01.wls.csv');
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/09.final_0.csv');
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/10.final_0_with_WLS.csv');
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/11.final_0_again.csv');

%%% nps
% ac_data = readtable('~/LOGS/swing/nps/test_04.csv', 'CommentStyle', '#'); meta = read_meta_swing("~/LOGS/swing/nps/test_04.csv");
ac_data = readtable('~/LOGS/swing/nps/test_01.csv', 'CommentStyle', '#'); meta = read_meta_swing("~/LOGS/swing/nps/test_01.csv");

% remove duplicates
[~, keep_idx] = unique(ac_data.timestamp, 'stable');
ac_data = ac_data(keep_idx, :);

% remove last line
ac_data(end, :) = [];
ac_data.timestamp = ac_data.timestamp - ac_data.timestamp(1);

Ts = median(diff(ac_data.timestamp)); 

disp("meta:"); disp(meta);
disp("Ts:"); disp(Ts);

%%
figure('Name', 'Translational');
plot_translational(ac_data);

%%
figure('Name', 'Eulers ZXY order');
plot_eulers(ac_data);

%%
figure('Name', 'Rates');
plot_rates(ac_data);

%%
figure('Name', 'Angular Acceleration');
plot_ang_accel(ac_data);

%%
figure('Name', 'Actuators');
plot_actuators(ac_data);

%%
figure('Name', 'Timings');
plot_timings(ac_data);

%%
figure('Name','Voltage');
plot(ac_data.timestamp, ac_data.voltage);

%% for old INDI
% figure('Name', 'Stabilization Loop');
% plot_stab_attitude(ac_data);