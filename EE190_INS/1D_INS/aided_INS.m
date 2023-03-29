clear; clc; close all;

load('SimplKF_data20170220A.mat');

%% Problem 5.1
%% Build a difference in time vector
dt = diff(t);

figure;
plot(dt, '.');
title("dt");
xlabel("samples");
ylabel("dt, seconds");
axis([0, inf, 7e-3, 9e-3])

figure;
plot(t, acc);
title("Acceleration");
xlabel("Time, t, (s)");
ylabel("x-acceleration, (m/s^2)");


%% Problem 5.2

%% Find an area of the data where the IMU is stationarty, average the data to remove the bias
b_acc = mean(acc(1:20 + 1)); %add one for base one indexing

%% Build the state vector for the system
x = zeros(3, size(t, 2));
x(3, 1:end) = b_acc;

%% Begin the INS
for ii = 2:size(dt, 2)
    %remove the bias off the measurment
    acc_unbiased = acc(ii) - x(3, ii);
    
    %integrate the acceleration to get velocity
    x(2, ii) = x(2, ii - 1) + acc_unbiased*dt(ii - 1);
    
    %integrate the velocity to get position
    x(1, ii) = x(1, ii - 1) + x(2, ii - 1)*dt(ii - 1) + 0.5*acc_unbiased*dt(ii - 1)^2;
end

x(1, end) = x(1, end - 1);
x(2, end) = x(2, end - 1);

figure;
subplot(2, 1, 1);
plot(t, x(1, :))
title("INS position, 5.2");
xlabel("time(s)");
ylabel("position (m)");
subplot(2, 1, 2); 
plot(t, x(2, :))
title("INS velocity, 5.2");
xlabel("time(s)");
ylabel("velocity (m/s)");

%% Problem 5.4

%% Get the number of samples in 20s worth of data and generate new time vector
sigma_n = 1e-3; %m/s/sqrt(s)
sigma_w = 1e-5; %m/s^2/sqrt(s)
tau = mean(dt);
f_s = (1 / mean(dt)); % IMU sampling rate
numPts = f_s * 20 + 1;
t_interval = t(1:numPts); 

% get the number of pts in 2s
numPts_2s = f_s * 2 + 1;
numPts_1s = f_s + 1;

clear x
P = zeros(3,3);
x = zeros(3, numPts);

%setup covariance propagation
phi = [1, tau, 0.5*tau^2;
       0, 1,   tau;
       0, 0,   1];

gamma = [0.5*tau^2, 0;
         tau, 0;
         0,   1];
Q = blkdiag(sigma_n^2, sigma_w^2);

%vector to extracts variance propagated
sigma_p = zeros(1, numPts_1s);
sigma_v = zeros(1, numPts_1s);

scale = 0;
for ii = 1:numPts_1s:numPts - numPts_1s
    %% calculate the bias for the IMU
    bias_acc = -mean(acc(ii - scale: ii + numPts_1s - 2 - scale));
    bias_var = var(acc(ii -scale : ii + numPts_1s - 2 - scale)) / f_s;
    P(3, 3) = bias_var;

    x(:, ii + numPts_1s - 1) = [0, 0, bias_acc]';
    for jj = ii + numPts_1s -1: ii + numPts_2s
        %remove the bias off the measurment
        x(3, jj) = bias_acc;
        acc_unbiased = acc(jj) + bias_acc;
    
        %integrate the acceleration to get velocity
        x(2, jj+1) = x(2, jj) + acc_unbiased*(1/f_s);
    
        %integrate the velocity to get position
        x(1, jj+1) = x(1, jj) + x(2, jj)*(1/f_s) + 0.5*acc_unbiased*(1/f_s)^2;

        if(ii == 1)
            P = phi*P*phi' + gamma*Q*gamma';
            sigma_p(jj - numPts_1s + ii) = sqrt(P(1,1));
            sigma_v(jj - numPts_1s + ii) = sqrt(P(2,2));
        end
    end
    scale = scale + 1;
end

x(1, end) = x(1, end - 1);
x(2, end) = x(2, end - 1);

figure;
subplot(4, 1, 1);
plot(t_interval(numPts_1s:end), x(1, numPts_1s:2501), '.')
title("INS position, 5.4");
xlabel("time(s)");
ylabel("position (m)");
axis([1, 20, -inf, inf])
subplot(4, 1, 2); 
plot(t_interval(numPts_1s:end), x(2, numPts_1s:2501), '.')
title("INS velocity, 5.4");
xlabel("time(s)");
ylabel("velocity (m/s)");
axis([1, 20, -inf, inf])
subplot(4, 1, 3); 
plot(t_interval(numPts_1s:end), x(3, numPts_1s:2501), '.')
title("INS bias, 5.4");
xlabel("time(s)");
ylabel("bias (m/s^2)");
axis([1, 20, -inf, inf])
subplot(4, 1, 4)
plot(t_interval(numPts_1s:end), acc(numPts_1s:2501), '.')
title("INS acc, 5.4");
xlabel("time(s)");
ylabel("acc (m/s^2)");
axis([1, 20, -inf, inf])

figure;
subplot(2,1,1)
for ii = numPts_1s:numPts_1s:numPts -numPts_1s
    t_bar = t_interval(ii + 1:ii + numPts_1s) - (((ii/f_s) -1) + (numPts_1s/f_s));
    plot(t_bar, x(1, ii: ii + numPts_1s -1), 'b');
    hold on;
end
plot(t_bar, sigma_p(1:end-1), 'k');
hold on
plot(t_bar, -sigma_p(1:end-1), 'k');
title("INS Error Position")
ylabel("Position Error Examples (m)")
xlabel("Integration Time t (s)")
axis([0, inf, -inf, inf])

subplot(2,1,2)
for ii = numPts_1s:numPts_1s:numPts -numPts_1s
    t_bar = t_interval(ii + 1:ii + numPts_1s) - (((ii/f_s) -1) + (numPts_1s/f_s));
    plot(t_bar, x(2, ii: ii + numPts_1s -1), 'b');
    hold on;
end
plot(t_bar, sigma_v(1:end-1), 'k');
hold on
plot(t_bar, -sigma_v(1:end-1), 'k');
title("INS Error Velocity")
xlabel("Integration Time t (s)")
ylabel("Velocity error Examples (m/s)");
axis([0, inf, -inf, inf])

clear x;
clear P;
clear sigma_v;
clear sigma_p;
close all;
%% Problem 5.5

sigma_n = 0.01;%m
H = [1, 0, 0];

%% Find an area of the data where the IMU is stationarty, average the data to remove the bias
b_acc = -mean(acc(1:20*f_s)); %add one for base one indexing

%% Build the state vector for the system
x_hat = zeros(3, size(t, 2));
P = zeros(3,3);
r = inf(1, size(t, 2));

sigma_p = zeros(1, size(t, 2));
sigma_v = zeros(1, size(t, 2));
sigma_b = zeros(1, size(t, 2));
sigma_r = inf(1, size(t, 2));


phi = [1, tau, 0.5*tau^2;
       0, 1,   tau;
       0, 0,   1];

gamma = [0.5*tau^2, 0;
         tau, 0;
         0,   1];

Q = blkdiag(sigma_n^2, sigma_w^2);

bias_var = var(acc(1:20 + 1)) / f_s;
P(3, 3) = bias_var;

x_hat(3, 1) = b_acc;

%remove in future
% close all;
%% Begin the INS
for ii = 2:size(dt, 2)
    if(mod(ii, (20*f_s)) ~= 0)
        x_hat(3, ii) = x_hat(3, ii - 1);

        %remove the bias off the measurment
        acc_unbiased = acc(ii) + x_hat(3, ii);
        
        %integrate the acceleration to get velocity
        x_hat(2, ii) = x_hat(2, ii - 1) + acc_unbiased*dt(ii - 1);
        
        %integrate the velocity to get position
        x_hat(1, ii) = x_hat(1, ii - 1) + x_hat(2, ii - 1)*dt(ii - 1) + 0.5*acc_unbiased*dt(ii - 1)^2;
    
        %Propagate the covariance
        P = phi*P*phi' + gamma*Q*gamma';
    else
        %prediction step for the kalman filter
        y_hat = H * x_hat(:, ii - 1);

        %compute the measurement residual, we know that at 20dt the IMU is
        %stationary
        r(ii) = 0 - y_hat;

        %calculate the covariance of the residual
        S = H*P*H' + sigma_n^2;

        %Kalman measurment update
        K = (P*H')/(S);
        P = P - K*(H*P);
        x_hat(:, ii) = x_hat(:, ii - 1) + K*r(ii);

        sigma_r(ii) = sqrt(S);
    end
    sigma_p(ii) = sqrt(P(1,1));
    sigma_v(ii) = sqrt(P(2,2));
    sigma_b(ii) = sqrt(P(3,3));
end

x_hat(1, end) = x_hat(1, end - 1);
x_hat(2, end) = x_hat(2, end - 1);

figure;
plot(t, r, '.');
hold on
plot(t, sigma_r, 'k.');
hold on
plot(t, -sigma_r, 'k.');
title("Residual");
xlabel("time (s)");
ylabel("Residual (m)")

figure;
subplot(3, 1, 1);
plot(t(1:12500), x_hat(1, 1:12500), '.')
hold on
plot(t(1:12500), x_hat(1, 1:12500) + sigma_p(1:12500), 'k.');
hold on
plot(t(1:12500), x_hat(1, 1:12500)-sigma_p(1:12500), 'k.')
title("Aided INS position, 5.5");
xlabel("time(s)");
ylabel("position (m)");
subplot(3, 1, 2); 
plot(t(1:12500), x_hat(2, 1:12500), '.')
hold on
plot(t(1:12500), x_hat(2, 1:12500)+sigma_v(1:12500), 'k.');
hold on
plot(t(1:12500), x_hat(2, 1:12500)-sigma_v(1:12500), 'k.')
title("Aided INS velocity, 5.5");
xlabel("time(s)");
ylabel("velocity (m/s)");
subplot(3,1,3)
plot(t, x_hat(3,:), '.');
hold on
plot(t, x_hat(3,:) + sigma_b, 'k.');
hold on
plot(t, x_hat(3,:)-sigma_b,'k.')
title("Aided INS bias, 5.5")
xlabel("time(s)");
ylabel("bias (m/s^2)");

P_p = (sigma_n^2/3)*t.^3 + (sigma_b(end)^2/4)*t.^4 + (sigma_w^2/20)*t.^5;
P_v = (sigma_n^2).*t + (sigma_b(end)^2).*t.^2 + (sigma_w^2/3).*t.^3;
figure;
subplot(2,1,1);
plot(t, P_p);
title("Position Error Covariance");
xlabel("time (s)");
ylabel("Error Covariance (m^2)");
subplot(2,1,2);
plot(t, P_v);
title("Velocity Error Covariance");
xlabel("time (s)");
ylabel("Error Covariance (m^2/s^2)");