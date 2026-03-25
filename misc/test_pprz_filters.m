% assumes ac_data loaded

fs = 1/median(diff(ac_data.timestamp));
t  = (ac_data.timestamp(1):1/fs:ac_data.timestamp(end))';

act_cmd = interp1(ac_data.timestamp(:), ac_data.act_cmd_1(:), t, 'linear');
act_state = interp1(ac_data.timestamp(:), ac_data.act_state_1(:), t, 'linear');

wc = 11;
a = exp(-wc/fs);

act_cmd_f = zeros(size(act_cmd));
act_cmd_f(1) = 0;

for k = 2:length(act_cmd)
    act_cmd_f(k) = a*act_cmd_f(k-1) + (1-a)*act_cmd(k);
end

% Plot comparison
figure; hold on; grid on;
plot(t, act_state, 'LineWidth',1,'DisplayName','act\_state');
plot(t, act_cmd_f, 'LineWidth',1,'DisplayName','act\_cmd\_f');
xlabel('Time [s]');
title('Actuator 1: measured vs expected first-order response');
legend;