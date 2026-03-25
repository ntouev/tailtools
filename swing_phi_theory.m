%% 
clear; 
% close all;
addpath('helpers/');
addpath('swing_plotters/');

%% log specifics and fixes
% ac_data = readtable('~/LOGS/swing/20260219_tuning/rpy.csv');
% tranges = [26 42]; % roll
% tranges = [38 60]; % pitch
% tranges = [64 76]; % thrust

% ac_data = readtable('~/LOGS/swing/20260304_stab/02.Kq2_Komega10.csv');
% tranges = [10 60];
% tranges = [15 35]; % pitch
% tranges = [39.5 45.5]; % roll
% ac_data = readtable('~/LOGS/swing/20260304_stab/08.Kq2.5Komega12.csv');
% tranges = [10 100];

% ac_data = readtable('~/LOGS/swing/20260306_more_stab/01.rigid_structure.csv');
% tranges = [10 120];
% ac_data = readtable('~/LOGS/swing/20260306_more_stab/03.training_dataset.csv');
% tranges = [10 90]; % roll picth
% tranges = [8 112]; % yaw

% ac_data = readtable('~/LOGS/swing/20260309_again_stab/01.wls.csv');
% tranges = [5 90];
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/09.final_0.csv');
% tranges = [5 55; 60 70];
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/10.final_0_with_WLS.csv');
% tranges = [3 38];
% ac_data = readtable('~/LOGS/swing/20260309_again_stab/11.final_0_again.csv');
% tranges = [5 66; 76 90];

% ac_data = readtable('~/LOGS/swing/20260310_optitrack/01_refinements_check.csv');
% tranges = [5 37];
% ac_data = readtable('~/LOGS/swing/20260310_optitrack/03.training.csv');
% tranges = [5 75];
% tranges = [5 25; 45 75]; % pitch
% tranges = [25 45]; % roll

% ac_data = readtable('~/LOGS/swing/20260311_training_data/03.training.csv');
% tranges = [15 120];
% ac_data = readtable('~/LOGS/swing/20260311_training_data/05.pitch90.csv');
% tranges = [6 56; 64 90];

% ac_data = readtable('~/LOGS/swing/20260312_vbatt/07.training.csv');
% tranges = [16 57; 61 130; 138 164];
% ac_data = readtable('~/LOGS/swing/20260312_vbatt/08.training_plus_vbattcode.csv');
% tranges = [8 72; 75 100; 103 122; 125 143; 145 151; 155 160];

ac_data = readtable('~/LOGS/swing/20260317_guidance/04.training_CT_3.5.csv');
tranges = [6 86];

% ac_data = readtable('~/LOGS/swing/20260319_quintic_again/14_yaw_training.csv');
% tranges = [35 50];

% remove duplicates
[~, keep_idx] = unique(ac_data.timestamp, 'stable');
ac_data = ac_data(keep_idx, :);

ac_data(end, :) = [];  % remove the last row
ac_data.timestamp = ac_data.timestamp - ac_data.timestamp(1);

%% t, datarange
fs = 100;
t = (0:1/fs:ac_data.timestamp(end))';

datarange = [];
for i = 1:size(tranges,1)
    trange = tranges(i,:);

    datarange_start = find(t > trange(1), 1, 'first') - 1;
    datarange_end = find(t > trange(2), 1, 'first') - 1;
    
    datarange = [datarange datarange_start:datarange_end];
end

%% interp (units: SI and pprz units)
qs = interp1(ac_data.timestamp, ac_data.qs, t, "linear", "extrap");
qx = interp1(ac_data.timestamp, ac_data.qx, t, "linear", "extrap");
qy = interp1(ac_data.timestamp, ac_data.qy, t, "linear", "extrap");
qz = interp1(ac_data.timestamp, ac_data.qz, t, "linear", "extrap");
quat = [qs qx qy qz];

p = interp1(ac_data.timestamp, ac_data.rates_p, t, "linear", "extrap");
q = interp1(ac_data.timestamp, ac_data.rates_q, t, "linear", "extrap");
r = interp1(ac_data.timestamp, ac_data.rates_r, t, "linear", "extrap");

v_n = interp1(ac_data.timestamp, ac_data.vel_n, t, "linear", "extrap");
v_e = interp1(ac_data.timestamp, ac_data.vel_e, t, "linear", "extrap");
v_d = interp1(ac_data.timestamp, ac_data.vel_d, t, "linear", "extrap");
vi = [v_n v_e v_d];
vb = quatrotate(quat, vi);

accx = interp1(ac_data.timestamp, ac_data.acc_x, t, "linear", "extrap");
accy = interp1(ac_data.timestamp, ac_data.acc_y, t, "linear", "extrap");
accz = interp1(ac_data.timestamp, ac_data.acc_z, t, "linear", "extrap");

act_cmd_TL = interp1(ac_data.timestamp, ac_data.act_cmd_1, t, "linear", "extrap");
act_cmd_TR = interp1(ac_data.timestamp, ac_data.act_cmd_2, t, "linear", "extrap");
act_cmd_BR = interp1(ac_data.timestamp, ac_data.act_cmd_3, t, "linear", "extrap");
act_cmd_BL = interp1(ac_data.timestamp, ac_data.act_cmd_4, t, "linear", "extrap");

vbatt = interp1(ac_data.timestamp, ac_data.voltage, t, "linear", "extrap");

%% 1st order actuator dynamics (pprz units)
G1 = tf(1, [1/11 1]);
actTL = lsim(G1, act_cmd_TL, t);
actTR = lsim(G1, act_cmd_TR, t);
actBR = lsim(G1, act_cmd_BR, t);
actBL = lsim(G1, act_cmd_BL, t);

%% filter with Butterworth
filter_freq = 4;
[b, a] = butter(2,filter_freq/(fs/2));

pf = filter(b, a, p, get_ic(b,a,p(1)));
qf = filter(b, a, q, get_ic(b,a,q(1)));
rf = filter(b, a, r, get_ic(b,a,r(1)));

vbf = filter(b, a, vb, get_ic(b,a,vb(1,:)));

accxf = filter(b, a, accx, get_ic(b,a,accx(1)));
accyf = filter(b, a, accy, get_ic(b,a,accy(1)));
acczf = filter(b, a, accz, get_ic(b,a,accz(1)));

actTLf = filter(b, a, actTL, get_ic(b,a,actTL(1)));
actTRf = filter(b, a, actTR, get_ic(b,a,actTR(1)));
actBRf = filter(b, a, actBR, get_ic(b,a,actBR(1)));
actBLf = filter(b, a, actBL, get_ic(b,a,actBL(1)));

vbattf = filter(b, a, vbatt, get_ic(b,a,vbatt(1)));

%% find derivatives of filtered values
pf_d = [zeros(1,1); diff(pf,1)]*fs;
qf_d = [zeros(1,1); diff(qf,1)]*fs;
rf_d = [zeros(1,1); diff(rf,1)]*fs;

pf_dd = [zeros(1,1); diff(pf_d,1)]*fs;
qf_dd = [zeros(1,1); diff(qf_d,1)]*fs;
rf_dd = [zeros(1,1); diff(rf_d,1)]*fs;

vbf_d = [zeros(1,3); diff(vbf,1)]*fs;

accxf_d = [zeros(1,1); diff(accxf,1)]*fs;
accyf_d = [zeros(1,1); diff(accyf,1)]*fs;
acczf_d = [zeros(1,1); diff(acczf,1)]*fs;

actTLf_d = [zeros(1,1); diff(actTLf,1)]*fs;
actTRf_d = [zeros(1,1); diff(actTRf,1)]*fs;
actBRf_d = [zeros(1,1); diff(actBRf,1)]*fs;
actBLf_d = [zeros(1,1); diff(actBLf,1)]*fs;

vbattf_d = [zeros(1,1); diff(vbattf,1)]*fs;

% weird derivatives
norm_vbf = sqrt(sum(vbf.*vbf, 2));
ddt_norm_vbf = (sum(vbf .* vbf_d, 2))./norm_vbf;

%% Fit Translational bx
Xtranx = [ddt_norm_vbf(datarange).*vbf(datarange,1) + norm_vbf(datarange).*vbf_d(datarange,1)]; 
Ytranx = accxf_d(datarange);

mdl_tranx = fitlm(Xtranx, Ytranx, "linear", 'Intercept', false);

%% Fit Translational by
Xtrany = [ddt_norm_vbf(datarange).*vbf(datarange,2) + norm_vbf(datarange).*vbf_d(datarange,2)]; 
Ytrany = accyf_d(datarange);

mdl_trany = fitlm(Xtrany, Ytrany, "linear", 'Intercept', false);

%% Fit Translational bz
% Xtranz = [ 2*actTLf(datarange).*actTLf_d(datarange), ...
%            2*actTRf(datarange).*actTRf_d(datarange), ...
%            2*actBRf(datarange).*actBRf_d(datarange), ... 
%            2*actBLf(datarange).*actBLf_d(datarange)]; 
% Ytranz = acczf_d(datarange);

Xtranz = [ddt_norm_vbf(datarange).*vbf(datarange,3) + norm_vbf(datarange).*vbf_d(datarange,3), ...
           vbattf(datarange).^2 .* 2.*actTLf(datarange).*actTLf_d(datarange), ...
           vbattf(datarange).^2 .* 2.*actTRf(datarange).*actTRf_d(datarange), ...
           vbattf(datarange).^2 .* 2.*actBRf(datarange).*actBRf_d(datarange), ... 
           vbattf(datarange).^2 .* 2.*actBLf(datarange).*actBLf_d(datarange)]; 
Ytranz = acczf_d(datarange);

% Xtranz = [ddt_norm_vbf(datarange).*vbf(datarange,3) + norm_vbf(datarange).*vbf_d(datarange,3), ...
%            2*actTLf(datarange).*actTLf_d(datarange), ...
%            2*actTRf(datarange).*actTRf_d(datarange), ...
%            2*actBRf(datarange).*actBRf_d(datarange), ... 
%            2*actBLf(datarange).*actBLf_d(datarange)]; 
% Ytranz = acczf_d(datarange);

mdl_tranz = fitlm(Xtranz, Ytranz, "linear", 'Intercept', false);

%% Fit Angular bx
% Xangx = [2*actTLf(datarange).*actTLf_d(datarange), ...
%          2*actTRf(datarange).*actTRf_d(datarange), ...
%          2*actBRf(datarange).*actBRf_d(datarange), ... 
%          2*actBLf(datarange).*actBLf_d(datarange)]; 
% Yangx = pf_dd(datarange);

Xangx = [vbattf(datarange).^2 .* 2.*actTLf(datarange).*actTLf_d(datarange), ...
         vbattf(datarange).^2 .* 2.*actTRf(datarange).*actTRf_d(datarange), ...
         vbattf(datarange).^2 .* 2.*actBRf(datarange).*actBRf_d(datarange), ... 
         vbattf(datarange).^2 .* 2.*actBLf(datarange).*actBLf_d(datarange)]; 
Yangx = pf_dd(datarange);

mdl_angx = fitlm(Xangx, Yangx, "linear", 'Intercept', false);

%% Fit Angular by
% Xangy = [2*actTLf(datarange).*actTLf_d(datarange), ...
%          2*actTRf(datarange).*actTRf_d(datarange), ...
%          2*actBRf(datarange).*actBRf_d(datarange), ... 
%          2*actBLf(datarange).*actBLf_d(datarange)]; 
% Yangy = qf_dd(datarange);

% Xangy = [ddt_norm_vbf(datarange).*vbf(datarange,1) + norm_vbf(datarange).*vbf_d(datarange,1), ...
%          2*actTLf(datarange).*actTLf_d(datarange), ...
%          2*actTRf(datarange).*actTRf_d(datarange), ...
%          2*actBRf(datarange).*actBRf_d(datarange), ... 
%          2*actBLf(datarange).*actBLf_d(datarange)]; 
% Yangy = qf_dd(datarange);

Xangy = [ddt_norm_vbf(datarange).*vbf(datarange,1) + norm_vbf(datarange).*vbf_d(datarange,1), ...
         vbattf(datarange).^2 .* 2.*actTLf(datarange).*actTLf_d(datarange), ...
         vbattf(datarange).^2 .* 2.*actTRf(datarange).*actTRf_d(datarange), ...
         vbattf(datarange).^2 .* 2.*actBRf(datarange).*actBRf_d(datarange), ... 
         vbattf(datarange).^2 .* 2.*actBLf(datarange).*actBLf_d(datarange)]; 
Yangy = qf_dd(datarange);

mdl_angy = fitlm(Xangy, Yangy, "linear", 'Intercept', false);

%% Fit Angular bz
% Xangz = [2*actTLf(datarange).*actTLf_d(datarange), ...
%          2*actTRf(datarange).*actTRf_d(datarange), ...
%          2*actBRf(datarange).*actBRf_d(datarange), ... 
%          2*actBLf(datarange).*actBLf_d(datarange)]; 
% Yangz = rf_dd(datarange);

Xangz = [vbattf(datarange).^2 .* 2.*actTLf(datarange).*actTLf_d(datarange), ...
         vbattf(datarange).^2 .* 2.*actTRf(datarange).*actTRf_d(datarange), ...
         vbattf(datarange).^2 .* 2.*actBRf(datarange).*actBRf_d(datarange), ... 
         vbattf(datarange).^2 .* 2.*actBLf(datarange).*actBLf_d(datarange)]; 
Yangz = rf_dd(datarange);

mdl_angz = fitlm(Xangz, Yangz, "linear", 'Intercept', false);

%% Plot
fprintf('Transl x: '); fprintf('%.3f ', mdl_tranx.Coefficients.Estimate); fprintf('\n');
fprintf('Transl y: '); fprintf('%.3f ', mdl_trany.Coefficients.Estimate); fprintf('\n');
fprintf('Transl z: '); fprintf('%.3f ', mdl_tranz.Coefficients.Estimate(1)); fprintf('\n');
fprintf('Transl z (x 10^-8): '); fprintf('%.3f ', mdl_tranz.Coefficients.Estimate(2:end)*1e8); fprintf('\n');

fprintf('Ang x  (x 10^-8): '); fprintf('%.2f ', mdl_angx.Coefficients.Estimate*1e8); fprintf('\n');
fprintf('Ang y: '); fprintf('%.3f ', mdl_angy.Coefficients.Estimate(1)); fprintf('\n');
fprintf('Ang y  (x 10^-8): '); fprintf('%.1f ', mdl_angy.Coefficients.Estimate(2:end)*1e8); fprintf('\n');
fprintf('Ang z  (x 10^-8): '); fprintf('%.2f ', mdl_angz.Coefficients.Estimate*1e8); fprintf('\n');

% fprintf('Accel x: '); fprintf('%.3f ', mdl_tranx.Coefficients.Estimate); fprintf('\n');
% fprintf('Accel y: '); fprintf('%.3f ', mdl_trany.Coefficients.Estimate); fprintf('\n');
% fprintf('Accel z: '); fprintf('%.9f ', mdl_tranz.Coefficients.Estimate); fprintf('\n');
% 
% fprintf('Ang x: '); fprintf('%.9f ', mdl_angx.Coefficients.Estimate); fprintf('\n');
% fprintf('Ang y: '); fprintf('%.9f ', mdl_angy.Coefficients.Estimate); fprintf('\n');
% fprintf('Ang z: '); fprintf('%.9f ', mdl_angz.Coefficients.Estimate); fprintf('\n');

figure('Name','Phi theory fit');
tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Ytranx, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_tranx.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[m/s^3]');
title('Jerk x');
legend('show');
hold off;

ax2 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yangx, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_angx.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[rad/s^3]');
title('Angular jerk x');
legend('show');
hold off;


ax3 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Ytrany, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_trany.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[m/s^3]');
title('Jerk y');
legend('show');
hold off;

ax4 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yangy, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_angy.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[rad/s^3]');
title('Angular jerk y');
legend('show');
hold off;

ax5 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Ytranz, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_tranz.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[m/s^3]');
title('Jerk z');
legend('show');
hold off;

ax6 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yangz, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_angz.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[rad/s^3]');
title('Angular jerk z');
legend('show');
hold off;

linkaxes([ax1, ax2, ax3, ax4, ax5, ax6],'x');

%%
figure('Name', 'Actuators');
plot_actuators(ac_data);

figure('Name', 'Translational');
plot_translational(ac_data);

figure('Name', 'body z accel');
plot(ac_data.timestamp, ac_data.acc_z);