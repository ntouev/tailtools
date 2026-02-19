clear;

homeDir = getenv('HOME');
addpath(genpath(fullfile(homeDir,'tails_plotters/swing_plotters')));

%%% real flight data
% ac_data = readtable('~/LOGS/swing/20260210_G_tuning/r-p-y-thr_20_30.csv');
ac_data = readtable('~/LOGS/swing/20260219_tuning/rpy.csv');
% ac_data = readtable('~/LOGS/swing/20260219_tuning/yaw_issue.csv');

%%% nps
% ac_data = readtable('~/LOGS/swing/nps/20260217-121852.csv');

ac_data(end, :) = [];
ac_data.timestamp = ac_data.timestamp - ac_data.timestamp(1);

%%
figure('Name', 'Stabilization Loop');
plot_stab_attitude(ac_data);

%%
figure('Name', 'Actuators');
plot_actuators(ac_data);