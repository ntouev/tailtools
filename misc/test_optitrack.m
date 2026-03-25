pos = ac_data.pos_n;
vel = ac_data.vel_n;
t   = ac_data.timestamp;

% pick only points where pos changes (20 Hz samples)
idx = [true; diff(pos) ~= 0];

pos_ds = pos(idx);
t_ds   = t(idx);

vel_num = gradient(pos_ds, t_ds);

figure;
plot(t_ds, vel_num, 'r--'); hold on;
plot(t, vel, 'b');
legend('Numerical differentiation', 'Measured');
xlabel('Time [s]');
ylabel('Velocity');
grid on;