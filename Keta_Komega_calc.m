%% new INDI (flatness-based)
close all;

wc = 11;

Komega_list = [14];
Keta_list   = [2 2.5 3.5 4 4.5];

G = tf(1, [1/wc 1], 'InputDelay', 0.0);
I = tf(1,[1 0]);

Gomega = G*I;

for i = 1:numel(Komega_list)
    Komega = Komega_list(i);

    Geta = feedback(Komega*Gomega, 1) * I;

    figure; hold on;

    lgd = strings(1, numel(Keta_list));

    for j = 1:numel(Keta_list)
        Keta = Keta_list(j);

        T = feedback(Keta*Geta, 1);
        step(T);
        lgd(j) = sprintf('K\\eta = %g', Keta);
    end

    title(sprintf('Step response, K\\omega = %g', Komega));
    legend(lgd, 'Location', 'best');
    hold off;
end


%% for old INDI
% close all;
% 
% wc = 11; % actuator bandwidth
% 
% G = tf(1, [1/wc 1], 'InputDelay', 0.0);
% I = tf(1,[1 0]); % integrator
% 
% Gomega = G*I;
% Komega = 7;
% Geta = feedback(Komega*Gomega, 1)*I;
% Keta = 2.6;
% step(feedback(Keta*Geta,1));
% stepinfo(feedback(Keta*Geta,1));
% grid on;
% 
% REF_RATE_x = Komega;
% REF_ERR_x = 2*Keta*Komega;

% conversions for rltool
% G1 = ss(Gomega);
% G2 = ss(I);