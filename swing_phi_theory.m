%% 
clear; 
% close all;
addpath('helpers/');

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
ac_data = readtable('~/LOGS/swing/20260309_again_stab/11.final_0_again.csv');
tranges = [5 66; 76 90];

% remove duplicates
[~, keep_idx] = unique(ac_data.timestamp, 'stable');
ac_data = ac_data(keep_idx, :);

ac_data(end, :) = [];  % remove the last row
ac_data.timestamp = ac_data.timestamp - ac_data.timestamp(1);

%% t, datarange
fs = 500;
t = (0:1/fs:ac_data.timestamp(end))';

datarange = [];
for i = 1:size(tranges,1)
    trange = tranges(i,:);

    datarange_start = find(t > trange(1), 1, 'first') - 1;
    datarange_end = find(t > trange(2), 1, 'first') - 1;
    
    datarange = [datarange datarange_start:datarange_end];
end

%% interp (units: SI and pprz units)
p = interp1(ac_data.timestamp, ac_data.rates_p, t, "linear", "extrap");
q = interp1(ac_data.timestamp, ac_data.rates_q, t, "linear", "extrap");
r = interp1(ac_data.timestamp, ac_data.rates_r, t, "linear", "extrap");

accz = interp1(ac_data.timestamp, ac_data.acc_z, t, "linear", "extrap");

act_cmd_TL = interp1(ac_data.timestamp, ac_data.act_cmd_1, t, "linear", "extrap");
act_cmd_TR = interp1(ac_data.timestamp, ac_data.act_cmd_2, t, "linear", "extrap");
act_cmd_BR = interp1(ac_data.timestamp, ac_data.act_cmd_3, t, "linear", "extrap");
act_cmd_BL = interp1(ac_data.timestamp, ac_data.act_cmd_4, t, "linear", "extrap");

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

acczf = filter(b, a, accz, get_ic(b,a,accz(1)));

actTLf = filter(b, a, actTL, get_ic(b,a,actTL(1)));
actTRf = filter(b, a, actTR, get_ic(b,a,actTR(1)));
actBRf = filter(b, a, actBR, get_ic(b,a,actBR(1)));
actBLf = filter(b, a, actBL, get_ic(b,a,actBL(1)));

%% find derivatives of filtered values
pf_d = [zeros(1,1); diff(pf,1)]*fs;
qf_d = [zeros(1,1); diff(qf,1)]*fs;
rf_d = [zeros(1,1); diff(rf,1)]*fs;

pf_dd = [zeros(1,1); diff(pf_d,1)]*fs;
qf_dd = [zeros(1,1); diff(qf_d,1)]*fs;
rf_dd = [zeros(1,1); diff(rf_d,1)]*fs;

acczf_d = [zeros(1,1); diff(acczf,1)]*fs;

actTLf_d = [zeros(1,1); diff(actTLf,1)]*fs;
actTRf_d = [zeros(1,1); diff(actTRf,1)]*fs;
actBRf_d = [zeros(1,1); diff(actBRf,1)]*fs;
actBLf_d = [zeros(1,1); diff(actBLf,1)]*fs;

%% Fit Accel z
% Xaccelz = [actTLf(datarange).^2, ...
%            actTRf(datarange).^2, ...
%            actBRf(datarange).^2, ... 
%            actBLf(datarange).^2]; 
% Yaccelz = acczf(datarange);

Xaccelz = [2*actTLf(datarange).*actTLf_d(datarange), ...
           2*actTRf(datarange).*actTRf_d(datarange), ...
           2*actBRf(datarange).*actBRf_d(datarange), ... 
           2*actBLf(datarange).*actBLf_d(datarange)]; 
Yaccelz = acczf_d(datarange);

mdl_accelz = fitlm(Xaccelz, Yaccelz, "linear", 'Intercept', false);

%% Fit Ang Accel x
% Xangx = [actTLf(datarange).^2, ...
%          actTRf(datarange).^2, ...
%          actBRf(datarange).^2, ... 
%          actBLf(datarange).^2]; 
% Yangx = pf_d(datarange);

Xangx = [2*actTLf(datarange).*actTLf_d(datarange), ...
         2*actTRf(datarange).*actTRf_d(datarange), ...
         2*actBRf(datarange).*actBRf_d(datarange), ... 
         2*actBLf(datarange).*actBLf_d(datarange)]; 
Yangx = pf_dd(datarange);

mdl_angx = fitlm(Xangx, Yangx, "linear", 'Intercept', false);

%% Fit Ang Accel y
% Xangy = [actTLf(datarange).^2, ...
%          actTRf(datarange).^2, ...
%          actBRf(datarange).^2, ... 
%          actBLf(datarange).^2]; 
% Yangy = pf_d(datarange);

Xangy = [2*actTLf(datarange).*actTLf_d(datarange), ...
         2*actTRf(datarange).*actTRf_d(datarange), ...
         2*actBRf(datarange).*actBRf_d(datarange), ... 
         2*actBLf(datarange).*actBLf_d(datarange)]; 
Yangy = qf_dd(datarange);

mdl_angy = fitlm(Xangy, Yangy, "linear", 'Intercept', false);

%% Fit Ang Accel z
% Xangz = [actTLf(datarange).^2, ...
%          actTRf(datarange).^2, ...
%          actBRf(datarange).^2, ... 
%          actBLf(datarange).^2]; 
% Yangz = pf_d(datarange);

Xangz = [2*actTLf(datarange).*actTLf_d(datarange), ...
         2*actTRf(datarange).*actTRf_d(datarange), ...
         2*actBRf(datarange).*actBRf_d(datarange), ... 
         2*actBLf(datarange).*actBLf_d(datarange)]; 
Yangz = rf_dd(datarange);

mdl_angz = fitlm(Xangz, Yangz, "linear", 'Intercept', false);

%% Plot
fprintf('Accel z (x 10^-8): '); fprintf('%.2f ', mdl_accelz.Coefficients.Estimate*1e8); fprintf('\n');

fprintf('Ang x   (x 10^-8): '); fprintf('%.1f ', mdl_angx.Coefficients.Estimate*1e8); fprintf('\n');
fprintf('Ang y   (x 10^-8): '); fprintf('%.0f ', mdl_angy.Coefficients.Estimate*1e8); fprintf('\n');
fprintf('Ang z   (x 10^-8): '); fprintf('%.1f ', mdl_angz.Coefficients.Estimate*1e8); fprintf('\n');

figure('Name','Phi theory fit');
tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yaccelz, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_accelz.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[m/s^2]');
title('DUMMY');
legend('show');
hold off;

ax2 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yangx, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_angx.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[rad/s^2]');
title('Angular Accel x');
legend('show');
hold off;


ax3 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yaccelz, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_accelz.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[m/s^2]');
title('DUMMY');
legend('show');
hold off;

ax4 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yangy, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_angy.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[rad/s^2]');
title('Angular Accel y');
legend('show');
hold off;

ax5 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yaccelz, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_accelz.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[m/s^2]');
title('Accel z');
legend('show');
hold off;

ax6 = nexttile;
hold on; grid on; zoom on;
plot(t(datarange), Yangz, '.', MarkerEdgeColor='b', DisplayName="Real data", LineWidth=1.5);
plot(t(datarange), mdl_angz.Fitted, '.', MarkerEdgeColor='r', DisplayName="Interpolated data", LineWidth=1.5);
xlabel('t [sec]');
ylabel('[rad/s^2]');
title('Angular Accel z');
legend('show');
hold off;

linkaxes([ax1, ax2, ax3, ax4, ax5, ax6],'x');

%%
figure('Name', 'Actuators');
plot_actuators(ac_data);