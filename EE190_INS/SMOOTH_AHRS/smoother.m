clear; clc;

%% Load in the Data
T = readtable("/home/kevinurrutia/INS/EE190_INS/data/Epson_G370_20230407_030457.csv");
%create vectors of the data
delta_v = [T.delta_v_x, T.delta_v_y, T.delta_v_z]';
delta_th = [T.delta_th_x, T.delta_th_y, T.delta_th_z]';
delta_th = delta_th * (pi/180); %covert to radians
t = seconds(T.system_time - T.system_time(1));

%throw away first sample since it is an outlier
delta_v = delta_v(:, 2:end);
delta_th = delta_th(:, 2:end);
test_delta_th = zeros(3, size(delta_th, 2));
t = t(2:end);

dt = diff(t);
tau = mean(dt);

%% Constants
fsamp = 125;
reset_period = 20;
tot_samps = size(delta_v, 2);
reset_samples = reset_period * fsamp;
N = 3.6 *(pi/180) * (1/3600);
Srg = N^2;
B = 0.9 * (pi/180) * (1/3600);
Tb = 60 / 1.89;
Sbg = (2*B^2*log(2)) / (pi*(0.4365)^2*Tb);

%% Level the IMU in the period where we know it is at rest
[ini_phi, ini_theta] = level(delta_v(:, 1:fsamp), fsamp);
ini_psi = 0; % we have no way of knowing from the data the orientation, we know because of experiment
ini_rpy = [ini_phi, ini_theta, ini_psi];

%% initialize AHRS
rpy = inf(3, tot_samps);
r = inf(3, tot_samps);
b = inf(3, tot_samps);
a = zeros(6, tot_samps);
c = zeros(3, tot_samps);
Phi_tensor = zeros(6, 6, tot_samps);
H_tensor = zeros(3, 6, tot_samps);
C_t_b_o = t2b(ini_rpy);
b_o = mean(delta_th(:, 1:fsamp), 2) / tau;
rpy(:, 1) = extract_rpy(C_t_b_o);
b(:, 1) = b_o;

%initialize the state covariance matrix
P = zeros(6,6);
Phi = zeros(6,6);
bias_var = var(delta_th(:, 1:fsamp)/tau, 0, 2);
sigma_rho = diag([0.01, 0.01, 0.01]);
P(1:3, 1:3) = sigma_rho^2*eye(3);
P(4:6, 4:6) = blkdiag(bias_var(1), bias_var(2), bias_var(3));

%setup measurement noise matrix
Q_tensor = zeros(6, 6, tot_samps);

%initialize the measurement covariance matrix
sig_meas = 1e-4; %[rad] error in measurement
R = sig_meas^2 * eye(3);
R_tensor = zeros(3, 3, tot_samps);

%Measurement Matrix
H = zeros(3, 6);
H(1:3, 1:3) = eye(3);

%setup smoothing matrices
A_dyn = [];
A_meas = [];
b_dyn = [];
b_meas = [];

C_b_t_o = C_t_b_o';
C_b_t_old = C_b_t_o;
C_b_t_new = C_b_t_old;

%setup delta vector to end smoothing process
delta = inf(6, 1); 
thresh_rho = 0.01; %the error between the attitude should be less than 0.01 rad
while norm(delta(1:3), 1) > thresh_rho
    for ii = 1:reset_samples:tot_samps - reset_samples
        %propagate the AHRS until the next measurement
        [C_b_t_hat, b_hat, rpy_hat, Phi_accum, a(:, ii), Q_accum] = AHRS(C_b_t_old, b(:, ii), delta_th(:, ii:ii + reset_samples), tau, Sbg, Srg);

        %construct measurement residuals
        c(:, ii) = measDiff(delta_v(:, ii: ii + reset_samples), C_b_t_hat, fsamp);

        %save tensors
        Phi_tensor(:, :, ii) = Phi_accum;
        H_tensor(:, :, ii) = H;
        Q_tensor(:, :, ii) = Q_accum;
        R_tensor(:, :, ii) = R;

        %save for plotting
        rpy(:, ii:ii + reset_samples) = rpy_hat;
        b(:, ii:ii + reset_samples) = b_hat;
        r(:, ii) = c(:, ii); 

        C_b_t_old = C_b_t_hat;
    end
    %clean up tensors
    span = find(~any(c, 1));
    c(:, span) = [];
    a(:, span) = [];
    Phi_tensor(:, :, span) = [];
    H_tensor(:, :, span) = [];
    Q_tensor(:, :, span) = [];
    R_tensor(:, :, span) = [];
    %Construct Smoothing matrices
    [A, B] = smoothMat(c, a, H_tensor, Phi_tensor, Q_tensor, R_tensor);
    delta = computeDelta(A, B);
    break
end

%% plotting
figure(1); clf;
subplot(3, 1, 1);
plot(t, rpy(1, :), ".");
grid on
title("Roll, Pitch, Yaw in Local Tangent Frame");
ylabel("Roll [rad]");
xlabel("time [s]");
subplot(3, 1, 2);
plot(t, rpy(2, :), ".");
grid on
ylabel("Pitch [rad]");
xlabel("time [s]");
subplot(3, 1, 3);
plot(t, rpy(3, :), ".");
grid on
ylabel("Yaw [rad]");
xlabel("time [s]");

figure(2); clf;
subplot(3, 1, 1);
plot(t, b(1, :));
grid on
title("Gyro Bias Estimate");
ylabel("bg_x [rad/s]");
xlabel("time [s]");
subplot(3, 1, 2);
plot(t, b(2, :));
grid on
title("Gyro Bias Estimate");
ylabel("bg_y [rad/s]");
xlabel("time [s]");
subplot(3, 1, 3);
plot(t, b(3, :));
grid on
title("Gyro Bias Estimate");
ylabel("bg_z [rad/s]");
xlabel("time [s]");

figure(3); clf;
subplot(3, 1, 1);
plot(t, r(1, :), '.');
grid on
title("Measurement Residuals")
ylabel("Attitude Error X [rad]")
xlabel("time [s]");
subplot(3, 1, 2);
plot(t, r(2, :), '.');
grid on
ylabel("Attitude Error Y [rad]")
xlabel("time [s]");
subplot(3, 1, 3);
plot(t, r(3, :), '.');
grid on
ylabel("Attitude Error Z [rad]")
xlabel("time [s]");

%% Functions
function delta = computeDelta(A, b)
    R = chol(A'*A);
    y = A' * b \ R';
    delta = y / R;
end

function [A, b] = smoothMat(c, a, H, Phi, Q, R)
    N = size(c, 2);
    measDim = size(c, 1);
    stateDim = size(Phi, 1);
    l = measDim * N * 2;
    w = stateDim * N + N * measDim;
    b = zeros(w, 1);
    A = zeros(w, l);

    for ii = 0:N-1
        sqrt_Q = sqrtm(Q(:, :, ii + 1));
        sqrt_R = sqrtm(R(:, :, ii + 1));

        G = -eye(stateDim);
        A(ii * stateDim + 1: ii * stateDim + stateDim, ii * stateDim + stateDim + 1: ii * stateDim + stateDim*2) = -sqrt_Q*G;
        A(ii * stateDim + 1: ii * stateDim + stateDim, ii * stateDim + 1: ii * stateDim + stateDim) = -sqrt_Q * Phi(:, :, ii + 1);

        A(ii * measDim + stateDim*N + 1: ii * measDim + stateDim *N + measDim, ii * stateDim + 1: ii * stateDim + stateDim) = sqrt_R*H(:, :, ii + 1);

        b(ii*stateDim + 1: ii * stateDim + stateDim, 1) = sqrt_Q * a(:, ii + 1);
        b(ii * measDim + stateDim * N + 1: ii * measDim + stateDim * N + measDim, 1) = sqrt_R * c(:, ii + 1);
    end
end

function Q = calcINSNoise(tau, Srg, Sbgd, C_b_i)
    Q = zeros(6,6);
    Q(1:3, 1:3) = C_b_i*Srg*eye(3)*C_b_i'*tau + (1/3)*C_b_i*Sbgd*eye(3)*C_b_i'*tau^3;
    Q(1:3, 4:6) = -(1/2)*C_b_i*Sbgd*eye(3)*tau^2;
    Q(4:6, 1:3) = -(1/2)*C_b_i'*Sbgd*eye(3)*tau^2;
    Q(4:6, 4:6) = Sbgd*tau*eye(3);
end

function r = measDiff(delta_v, C_b_t_hat, fsamp)
    [phi, theta] = level(delta_v, fsamp);
    psi = 0;
    rpy_meas = [phi, theta, psi];
    C_b_t_o = t2b(rpy_meas)';
    deltaC = C_b_t_o * C_b_t_hat';
    rho_skew = deltaC - eye(3);
    r = inv_skew_symmetric(rho_skew);
end
function [C_b_t_hat, b_hat, rpy_hat, Phi_accum, a, Q_accum] = AHRS(C_b_t_hat_old, b, alpha, tau, Sbg, Srg)
    Phi_accum = eye(6);
    Q_accum = zeros(6);
    b_hat = zeros(3, size(alpha, 2));
    rpy_hat = zeros(3, size(alpha, 2));
    C_b_t_hat = C_b_t_hat_old;
    b_hat(:, 1) = b;
    rpy_hat(:, 1) = extract_rpy(C_b_t_hat);
    for jj = 2:size(alpha, 2)
        %propagate bias
        b_hat(:, jj) = b_hat(:, jj - 1);

        %remove the bias from the measurement
        alpha(:, jj) = alpha(:, jj) - b_hat(:, jj) * tau;

        %propagate the rotation matrix;
        C_b_t_hat = attitude_update(C_b_t_hat, alpha(:, jj));

        %extract rpy
        rpy_hat(:, jj) = extract_rpy(C_b_t_hat);

        %calculate state transition matrix
        Phi = calcStateTransition(C_b_t_hat, tau);
        Q = calcINSNoise(tau, Srg, Sbg, C_b_t_hat);

        Q_accum = Q + Q_accum;
        Phi_accum = Phi * Phi_accum;
    end
    
    %return the error between the old state and the propagated one
    delta_C = C_b_t_hat * C_b_t_hat_old';
    rho_skew = delta_C - eye(3);
    rho = inv_skew_symmetric(rho_skew);
    b_err = b - b_hat(:, end);
    a = [rho; b_err];
end

function  [phi, theta] = level(delta_v, fsamp)
    theta_arr = zeros(1, 1*fsamp);
    phi_arr = zeros(1, 1*fsamp);
    
    for ii = 1:round(1*fsamp)
        theta_arr(ii) = atan(delta_v(1, ii) / sqrt((delta_v(2, ii)^2  + (delta_v(3, ii))^2)));
        phi_arr(ii) = atan2(-delta_v(2, ii), -delta_v(3, ii));
    end

    phi = mean(phi_arr);
    theta = mean(theta_arr);
end

function C_i_b = t2b(delta_th)
   phi = delta_th(1);
   theta = delta_th(2);
   psi = delta_th(3);

   C_i_b(1, 1) = cos(theta)*cos(psi);
   C_i_b(1, 2) = cos(theta)*sin(psi);
   C_i_b(1, 3) = -sin(theta);
   C_i_b(2, 1) = -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi);
   C_i_b(2, 2) = cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi);
   C_i_b(2, 3) = sin(phi)*cos(theta);
   C_i_b(3, 1) = sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
   C_i_b(3, 2) = -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
   C_i_b(3, 3) = cos(phi)*cos(theta);
end

function rpy = extract_rpy(C_b_i) 
    phi = atan2(C_b_i(3,2), C_b_i(3,3));
    theta = -atan((C_b_i(3,1))/(sqrt(1 - (C_b_i(3,1))^2)));
    psi = atan2(C_b_i(2,1), C_b_i(1,1));


    rpy = [(phi), (theta), (psi)]';
end

function C_b_i_new = attitude_update(C_b_i_old, alpha)
    Omega_ib_b = Skew_symmetric(alpha);

    C_b_i_new = C_b_i_old*(eye(3) + Omega_ib_b);
end

function A = Skew_symmetric(a)
    A = [    0, -a(3),  a(2);...
      a(3),     0, -a(1);...
     -a(2),  a(1),     0];
end

function phi = calcStateTransition(C_b_i, tau)
    phi = zeros(6,6);

    phi(1:3, 1:3) = eye(3);
    phi(1:3, 4:6) = -C_b_i*tau;
    phi(4:6, 1:3) = zeros(3,3);
    phi(4:6, 4:6) = eye(3);

end

function x = inv_skew_symmetric(R)
    x = zeros(3, 1);

    x(1) = mean([-R(2, 3),R(3, 2)]);
    x(2) = mean([-R(3, 1), R(1, 3)]);
    x(3) = mean([-R(1, 2), R(2, 1)]);
end