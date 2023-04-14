clear; clc; close all;

%% load in the data
T = readtable("Epson_G370_20230407_030457.csv");
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

%% constants
fsamp = 125; %Hz
% fcut = fsamp/2; %Hz cutoff frequency
fcut = 52.4; %Hz
reset_period = 1; %s, IMU is placed back at rest at this time
reset_samples = reset_period * fsamp; %number of samples from last reset period
plot_until = 18; %s, number of seconds to plot at end of run
tot_samps = plot_until * fsamp;%size(delta_v, 2); % number of samples in the experiment
Srg = ((0.6*(1/60)*(pi/180))*sqrt(tau))^2; %degrees/sqrt(s), angular random walk
% Sbgd = ((0.8*(1/3600)*(pi/80))*3; %degree/hour, gyro bias
B = ((0.8*(1/3600)*(pi/80))^2*pi)/(2*log(2));
% Sbgd = (B^2)/(2*pi*fcut);
Sbgd = (2*B^2*log(2))/(pi*(0.4365)^2*(70/1.89));

%% Level the IMU in the period where we know it is at rest
[ini_phi, ini_theta] = level(delta_v(:, 1:fsamp), fsamp);
ini_psi = 0; % we have no way of knowing from the data the orientation, we know because of experiment
ini_rpy = [ini_phi, ini_theta, ini_psi];

%% initialize AHRS
rpy = zeros(3, tot_samps);
b = zeros(3, tot_samps);
C_t_b_o = t2b(ini_rpy); % initial rot matrix from body to NED should be identity since ini_rpy is known to be about [0, 0, 0]'
bias_th_hat = mean(delta_th(:,1:round(fsamp)),2);
rpy(:, 1) = extract_rpy(C_t_b_o');
b(:, 1) = bias_th_hat;

P = zeros(6,6);
Phi = zeros(6,6);
bias_var = std(delta_th(:,1:round(fsamp)), 0, 2);
sigma_rho = (pi/4);
sigma_b = 1e-5;
P(1:3, 1:3) = sigma_rho^2*eye(3); %A very large error since we are uncertain about the state
% P(4:6, 4:6) = blkdiag(bias_var(1)^2, bias_var(2)^2, bias_var(3)^2);
P(4:6, 4:6) = sigma_b^2*eye(3);

Q = zeros(6,6);

%use the following array to extract the covariance being propagated for
%the state
var_rho = zeros(3, tot_samps); %covariance for the attitude error
var_b = zeros(3, tot_samps); %covariance for the attitude increment bias

var_rho(:, 1) = [P(1,1), P(2,2), P(3,3)]';
var_b(:, 1) = [P(4,4), P(5,5), P(6,6)]';

%initialize residual vector
r = inf(3, tot_samps);
sig_meas = 1e-4; %[rad], error in model
R = sig_meas^2 * eye(3);

%use the following arrays to propagate the residual covariance
var_r = inf(3, tot_samps);

%Measurement Matrix
H = zeros(3, 6);
H(1:3, 1:3) = eye(3);

K_plot = zeros(18, tot_samps);
K = zeros(6,3);

C_b_t_o = C_t_b_o';
C_b_t_old = C_b_t_o;
C_b_t_new = C_b_t_old;
for ii = 2:tot_samps
    tau = t(ii) - t(ii - 1);
    if(mod(ii, (reset_period*fsamp)) ~= 0)
        %propagate bias
        b(:, ii) = b(:, ii - 1);
    
        %remove bias from measurement
        delta_th(:, ii) = delta_th(:, ii) - b(:, ii)*tau;
    
        %propagate rotation matrix using delta_th
        C_b_t_new = attitude_update(C_b_t_old, delta_th(:, ii));
    
        %extract rpy for plotting
        rpy(:, ii) = extract_rpy(C_b_t_new);
    
        %Calculate state transition matrix
        Phi = calcStateTransition(C_b_t_new, tau);
    
        %Calculate the noise covariance matrix
        Q = calcINSNoise(tau, Srg, Sbgd, C_b_t_new);
    
        P = Phi*P*Phi' + Q;

        %update rotation matrix
        C_b_t_old = C_b_t_new;
    else
        %re-level the IMU
        [phi, theta] = level(delta_v(:, ii:ii + 1*fsamp), fsamp);
        psi = 0; % we have no way of knowing from the data the orientation, we know because of experiment
        rpy_meas = [phi, theta, psi];
        C_b_t_o = t2b(rpy_meas)';
        %measurement prediction step
        delta_C = C_b_t_o' * C_b_t_old;
        rho_skew = delta_C - eye(3);
        rho = inv_skew_symmetric(rho_skew);

        r(:, ii) = rho;

        %residual covariance
        S = H*P*H' + R;

        %Kalman measurement update
        K = (P*H')/S;

        %covariance update
        P = P - K*(H*P);

        %state update
        error_state = K*r(:,ii);

        %update rotation matrix
        C_b_t_new = C_b_t_old*(eye(3) + Skew_symmetric(error_state(1:3)))';
        C_b_t_old = C_b_t_new;

        %extract rpy for plotting
        rpy(:, ii) = extract_rpy(C_b_t_old);

        b(:, ii) = b(:, ii - 1) + error_state(4:6);

        var_r(:, ii) = [S(1,1), S(2,2), S(3,3)]';

        %store K for plotting
        K_plot(:, ii) = reshape(K, 1, []);
        
    end
    %store K for plotting
    K_plot(:, ii) = reshape(K, 1, []);
    var_rho(:, ii) = [P(1,1), P(2,2), P(3,3)]';
    var_b(:, ii) = [P(4,4), P(5,5), P(6,6)]';
end

axis_plot = 0.009;

%% Plot Results
figure;
subplot(3,1,1);
plot(t(1:plot_until*fsamp),real(rpy(1, 1:plot_until*fsamp)), '.');
hold on 
plot(t(1:plot_until*fsamp),real(rpy(1, 1:plot_until*fsamp)) - sqrt(var_rho(1, 1:plot_until*fsamp)), 'k.');
hold on 
plot(t(1:plot_until*fsamp),real(rpy(1, 1:plot_until*fsamp)) + sqrt(var_rho(1, 1:plot_until*fsamp)), 'k.');
title("Roll (rad)");
xlabel("time (s)");
ylabel("Roll (rad)");
axis([0, inf, -axis_plot, axis_plot])
subplot(3,1,2);
plot(t(1:plot_until*fsamp),real(rpy(2, 1:plot_until*fsamp)), '.');
hold on 
plot(t(1:plot_until*fsamp),real(rpy(2, 1:plot_until*fsamp)) - sqrt(var_rho(2, 1:plot_until*fsamp)), 'k.');
hold on 
plot(t(1:plot_until*fsamp),real(rpy(2, 1:plot_until*fsamp)) + sqrt(var_rho(2, 1:plot_until*fsamp)), 'k.');
title("Pitch (rad)");
xlabel("time (s)");
ylabel("Pitch (rad)");
axis([0, inf, -axis_plot, axis_plot])
subplot(3,1,3);
plot(t(1:plot_until*fsamp),real(rpy(3, 1:plot_until*fsamp)), '.');
hold on 
plot(t(1:plot_until*fsamp),real(rpy(3, 1:plot_until*fsamp)) - sqrt(var_rho(3, 1:plot_until*fsamp)), 'k.');
hold on 
plot(t(1:plot_until*fsamp),real(rpy(3, 1:plot_until*fsamp)) + sqrt(var_rho(3, 1:plot_until*fsamp)), 'k.');
title("Yaw (rad)");
xlabel("t (s)");
ylabel("Yaw (rad)");
axis([0, inf, -axis_plot, axis_plot])

figure;
subplot(3,1,1);
plot(t(1:plot_until*fsamp), b(1, 1:plot_until*fsamp), '.');
hold on
plot(t(1:plot_until*fsamp), b(1, 1:plot_until*fsamp) - sqrt(var_b(1, 1:plot_until*fsamp)), '.k');
hold on
plot(t(1:plot_until*fsamp), b(1, 1:plot_until*fsamp) + sqrt(var_b(1, 1:plot_until*fsamp)), '.k');
title("b g_x (m)");
xlabel("time (s)");
ylabel("b g_x (m)");
subplot(3,1,2);
plot(t(1:plot_until*fsamp), b(2, 1:plot_until*fsamp), '.');
hold on
plot(t(1:plot_until*fsamp), b(2, 1:plot_until*fsamp) - sqrt(var_b(2, 1:plot_until*fsamp)), '.k');
hold on
plot(t(1:plot_until*fsamp), b(2, 1:plot_until*fsamp) + sqrt(var_b(2, 1:plot_until*fsamp)), '.k');
title("b g_y (m)");
xlabel("time (s)");
ylabel("b g_y (m)");
subplot(3,1,3);
plot(t(1:plot_until*fsamp), b(3, 1:plot_until*fsamp), '.');
hold on
plot(t(1:plot_until*fsamp), b(3, 1:plot_until*fsamp) - sqrt(var_b(3, 1:plot_until*fsamp)), '.k');
hold on
plot(t(1:plot_until*fsamp), b(3, 1:plot_until*fsamp) + sqrt(var_b(3, 1:plot_until*fsamp)), '.k');
title("b g_z (m)");
xlabel("t (s)");
ylabel("b g_z (m)");

axis_plot = 0.001;
figure;
%plot residual vs covariance of residual 
subplot(3, 1, 1);
plot(t(1:tot_samps), r(1, 1:tot_samps), 'b.');
hold on
plot(t(1:tot_samps), -sqrt(var_r(1, 1:tot_samps)), 'k.');
hold on
plot(t(1:tot_samps), sqrt(var_r(1, 1:tot_samps)), 'k.');
title("Error Attitude X");
xlabel("time");
ylabel("Error Attitude X");
axis([0, inf, -axis_plot, axis_plot])
subplot(3, 1, 2);
plot(t(1:tot_samps), r(2, 1:tot_samps), 'b.');
hold on
plot(t(1:tot_samps), -sqrt(var_r(2, 1:tot_samps)), 'k.');
hold on
plot(t(1:tot_samps), sqrt(var_r(2, 1:tot_samps)), 'k.');
title("Error Attitude Y");
xlabel("time");
ylabel("Error Attitude Y");
axis([0, inf, -axis_plot, axis_plot])
subplot(3, 1, 3);
plot(t(1:tot_samps), r(3, 1:tot_samps), 'b.');
hold on
plot(t(1:tot_samps), -sqrt(var_r(3, 1:tot_samps)), 'k.');
hold on
plot(t(1:tot_samps), sqrt(var_r(3, 1:tot_samps)), 'k.');
title("Error Attitude Z");
xlabel("time");
ylabel("Error Attitude Z");
axis([0, inf, -axis_plot, axis_plot])

% figure;
% %plot all of the elements of K
% subplot(3,1,1);
% plot(t(1:tot_samps), K_plot(1, 1:tot_samps));
% title("K_1");
% subplot(3,1,2);
% plot(t(1:tot_samps), K_plot(2, 1:tot_samps));
% title("K_2");
% subplot(3,1,3);
% plot(t(1:tot_samps), K_plot(2, 1:tot_samps));
% title("K_3");
% figure;
% %plot all of the elements of K
% subplot(3,1,1);
% plot(t(1:tot_samps), K_plot(4, 1:tot_samps));
% title("K_4")
% subplot(3,1,2);
% plot(t(1:tot_samps), K_plot(5, 1:tot_samps));
% title("K_5")
% subplot(3,1,3);
% plot(t(1:tot_samps), K_plot(6, 1:tot_samps));
% title("K_6")
% figure;
% %plot all of the elements of K
% subplot(3,1,1);
% plot(t(1:tot_samps), K_plot(7, 1:tot_samps));
% title("K_7")
% subplot(3,1,2);
% plot(t(1:tot_samps), K_plot(8, 1:tot_samps));
% title("K_8")
% subplot(3,1,3);
% plot(t(1:tot_samps), K_plot(9, 1:tot_samps));
% title("K_9")
% figure;
% %plot all of the elements of K
% subplot(3,1,1);
% plot(t(1:tot_samps), K_plot(10, 1:tot_samps));
% title("K_{10}")
% subplot(3,1,2);
% plot(t(1:tot_samps), K_plot(11, 1:tot_samps));
% title("K_{11}")
% subplot(3,1,3);
% plot(t(1:tot_samps), K_plot(12, 1:tot_samps));
% title("K_{12}")
% figure;
% %plot all of the elements of K
% subplot(3,1,1);
% plot(t(1:tot_samps), K_plot(13, 1:tot_samps));
% title("K_{13}")
% subplot(3,1,2);
% plot(t(1:tot_samps), K_plot(14, 1:tot_samps));
% title("K_{14}")
% subplot(3,1,3);
% plot(t(1:tot_samps), K_plot(15, 1:tot_samps));
% title("K_{15}")
% figure;
% %plot all of the elements of K
% subplot(3,1,1);
% plot(t(1:tot_samps), K_plot(16, 1:tot_samps));
% title("K_{16}")
% subplot(3,1,2);
% plot(t(1:tot_samps), K_plot(17, 1:tot_samps));
% title("K_{17}")
% subplot(3,1,3);
% plot(t(1:tot_samps), K_plot(18, 1:tot_samps));
% title("K_{18}")

%% Functions
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

function C_t_b = t2b(delta_th)
   phi = delta_th(1);
   theta = delta_th(2);
   psi = delta_th(3);

   C_t_b(1, 1) = cos(theta)*cos(psi);
   C_t_b(1, 2) = cos(theta)*sin(psi);
   C_t_b(1, 3) = -sin(theta);
   C_t_b(2, 1) = -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi);
   C_t_b(2, 2) = cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi);
   C_t_b(2, 3) = sin(phi)*cos(theta);
   C_t_b(3, 1) = sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi);
   C_t_b(3, 2) = -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi);
   C_t_b(3, 3) = cos(phi)*cos(theta);
end

function C_b_t_new = attitude_update(C_b_t_old, alpha)
    Omega_ib_b = Skew_symmetric(alpha);

    C_b_t_new = C_b_t_old*(eye(3) + Omega_ib_b);
end

function A = Skew_symmetric(a)
    A = [    0, -a(3),  a(2);...
      a(3),     0, -a(1);...
     -a(2),  a(1),     0];
end


function rpy = extract_rpy(C_b_t) 
    phi = atan2(C_b_t(3,2), C_b_t(3,3));
    theta = -atan((C_b_t(3,1))/(sqrt(1 - (C_b_t(3,1))^2)));
    psi = atan2(C_b_t(2,1), C_b_t(1,1));


    rpy = [real(phi), real(theta), real(psi)]';
end

function phi = calcStateTransition(C_b_t, tau)
    phi = zeros(6,6);

    phi(1:3, 1:3) = eye(3);
    phi(1:3, 4:6) = C_b_t*tau;
    phi(4:6, 1:3) = zeros(3,3);
    phi(4:6, 4:6) = eye(3);

end

function Q = calcINSNoise(tau, Srg, Sbgd, C_b_t)
    Q = zeros(6,6);
    Q(1:3, 1:3) = C_b_t*Srg*eye(3)*C_b_t'*tau + (1/3)*C_b_t*Sbgd*eye(3)*C_b_t'*tau^3;
    Q(1:3, 4:6) = (1/2)*C_b_t*Sbgd*eye(3)*tau^2;
    Q(4:6, 1:3) = (1/2)*C_b_t'*Sbgd*eye(3)*tau^2;
    Q(4:6, 4:6) = Sbgd*tau*eye(3);
end

function x = inv_skew_symmetric(R)
    x = zeros(3, 1);

    x(1) = mean([-R(2, 3),R(3, 2)]);
    x(2) = mean([-R(3, 1), R(1, 3)]);
    x(3) = mean([-R(1, 2), R(2, 1)]);
end
