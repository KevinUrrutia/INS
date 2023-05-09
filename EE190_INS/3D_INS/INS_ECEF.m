clear; clc;

%% load in the data
T = readtable("/home/kevinurrutia/INS/EE190_INS/data/Epson_G370_20230407_030457.csv");
%create vector from the data
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
reset_period = 1; %s, IMU is placed back at initial position and is at rest
reset_samples =  reset_period*fsamp;
% plot_until = size(delta_v, 2)*tau;
plot_until = 18*tau;
% tot_samps = size(delta_v, 2);
tot_samps = 18*fsamp;
lla_o = [33.979130, -117.372570, 827]';
r_e_o = LLA_to_ECEF(lla_o);
g_e_o = Gravity_ECEF(r_e_o);
Srg = ((1e-3)*(pi/180))^2;
B = ((0.9*(1/3600)*(pi/180)))/0.664; %(rad/s)/0.664 = B
Tb = 32; %s
Sbg = (2*B^2*log(2))/(pi*(0.4365)^2*(Tb));
G = 9.80065;%[m/s^2]
t_b = 100; %[s] flat region of ASD curve
Sra = ((50e-6)*G)^2; %N ([m/s^2]^2)/Hz
Tb = t_b / 1.89;
B = ((8e-6)*G)/0664;
Sba = (2*(B)^2*log(2)) / (pi*(0.4365)^2*(Tb));


%% Level IMU
[ini_phi, ini_theta] = level(delta_v(:, 1:fsamp), fsamp);
ini_psi = 0; % we have no way of knowing from the data the orientation, we know because of experiment
ini_rpy = [ini_phi, ini_theta, ini_psi];

%Initial rotation matrices
C_t_b_o = t2b(ini_rpy);
C_t_e_o = t2e(lla_o(1), lla_o(2));
C_b_e_o = C_t_e_o * C_t_b_o';

%% Calculate Biases
b_omega_hat = mean(delta_th(:, 1:round(fsamp)), 2)/tau;
b_f_hat = mean(delta_v(:, 1:round(fsamp)), 2)/tau + C_b_e_o'*g_e_o;

b_hat = zeros(6, tot_samps);
b_hat(1:3) = b_f_hat;
b_hat(4:6) = b_omega_hat;

%% Initialize State Covariance Matrix
P = zeros(15, 15);
b_omega_var = var(delta_th(:,1:round(fsamp))/tau, 0, 2);
b_f_var = var(delta_v(:,1:round(fsamp))/tau, 0, 2);
var_rho = diag([.1, .1, 0.01].^2); %[rad]
sigma_v_e = 0.1; %[m/s]
sigma_r_e = 0.1; %[m]
P(1:3, 1:3) = var_rho;
P(4:6, 4:6) = sigma_v_e^2*eye(3);
P(7:9, 7:9) = sigma_r_e^2*eye(3);
P(10:12, 10:12) = diag(b_f_var);
P(13:15, 13:15) = diag(b_omega_var);

%% Initialize Meausurement Noise Covariance Matrix
sigma_meas = 1e-4;
R = sigma_meas^2*eye(9);

%% Initialize Measurement Observation Matrix
H = zeros(9, 15);
H(1:9, 1:9) = eye(9);

%% Preallocate Save Space
rpy = zeros(3, tot_samps);
v_e = zeros(3, tot_samps);
r_e = zeros(3, tot_samps);
g_t = zeros(3, tot_samps);
resid = inf(9, tot_samps);
var_rho = zeros(3, tot_samps);
var_v_e = zeros(3, tot_samps);
var_r_e = zeros(3, tot_samps);
var_b = zeros(6, tot_samps);
var_resid = inf(9, tot_samps);

r_e(:, 1) = r_e_o;
rpy(:, 1) = extract_rpy(C_t_b_o');
g_t(:, 1) = C_t_e_o'*g_e_o;
var_rho(:, 1) = [P(1,1), P(2,2), P(3,3)]';
var_v_e(:, 1) = [P(4,4), P(5,5), P(6,6)]';
var_r_e(:, 1) = [P(7,7), P(8,8), P(9,9)]';
var_b(:, 1) = [P(10,10), P(11,11), P(12,12), P(13,13), P(14,14), P(15,15)]';

%% INS Propagation
C_b_e_old = C_b_e_o;
C_b_e_new = C_b_e_o;
g_e = g_e_o;
for ii = 2:tot_samps
    tau = t(ii) - t(ii - 1);

    %propagate bias
    b_hat(:, ii) = b_hat(:, ii - 1);

    %remove biases from measurements
    delta_v(:, ii) = delta_v(:, ii) - b_hat(1:3, ii)*tau;
    delta_th(:, ii) = delta_th(:, ii) - b_hat(4:6, ii)*tau + + C_b_e_old'*g_e*tau;

    %propagate the rotation matrix
    C_b_e_new = attitude_update(C_b_e_old, delta_th(:, ii), tau);

    %integrated specific force rotation
    nu_e = int_specific_force_rot(delta_v(:, ii), C_b_e_old, C_b_e_new);

    %velocity update
    v_e(:, ii) = velocity_update(v_e(:, ii - 1), nu_e, g_e, tau);

    %position update
    r_e(:, ii) = position_update(r_e(:, ii - 1), v_e(:, ii - 1), v_e(:, ii), tau);

    %update gravity
    g_e = Gravity_ECEF(r_e(:, ii));

    %calculate lat, lon, alt
    lla = ECEF_to_LLA(r_e(:, ii));

    %calculate Phi
    Phi = calcStateTransition(C_b_e_new, delta_v(:, ii), r_e(:, ii), g_e, lla(1), tau);

    %Calculate Process Noise covariance 
    Q_d = calcINSNoise(tau, Sra, Srg, Sba, Sbg);

    %Propagate State Covariance Matrix
    P = Phi*P*Phi' + Q_d;

    %Earth to tangent rotation matrix
    C_t_e = t2e(lla(1), lla(2));

    %update rotation matrix
    C_b_e_old = C_b_e_new;

    %% Mesurement Update
    if(mod(ii, (reset_period*fsamp)) == 0)
        %re-level the IMU
        [phi, theta] = level(delta_v(:, ii - fsamp+1:ii), fsamp);
        psi = 0;
        rpy_meas = [phi, theta, psi];
        C_b_t = t2b(rpy_meas)';
        C_b_e_o = C_t_e*C_b_t;
        delta_C = C_b_e_o * C_b_e_old;
        rho_skew = delta_C - eye(3);
        rho = inv_skew_symmetric(rho_skew);

        %error in velcoity
        error_v = -v_e(:, ii);

        %error in position
        error_r = r_e_o - r_e(:, ii);

        %residual constructrion
        resid(:, ii) = [rho; error_v; error_r];

        %residual covariance
        S = H*P*H' + R;

        %Kalman Gain
        K = (P*H')/S;

        %covariance update
        P = P - K*(H*P);

        %state update
        error_state = K*resid(:, ii);

        %update rotation matrix
        Pc        = Skew_symmetric(error_state(1:3));
        C_corr    = (eye(3) + Pc);
        C_b_e_new = C_corr * C_b_e_old;    
        C_b_e_old = C_b_e_new; 

        %update velocity
        v_e(:, ii) = v_e(:, ii) + error_state(4:6);

        %update position
        r_e(:, ii) = r_e(:, ii) + error_state(7:9);

        %update bias
        b_hat(:, ii) = b_hat(:, ii) + error_state(10:15);

        %update gravity
        g_e = Gravity_ECEF(r_e(:, ii));
    
        %calculate lat, lon, alt
        lla = ECEF_to_LLA(r_e(:, ii));

        %Earth to tangent rotation matrix
        C_t_e = t2e(lla(1), lla(2));

        %Extract residual covariance for plotting
        var_resid(:, ii) = [S(1,1), S(2,2), S(3,3), S(4,4), S(5,5), S(6,6), S(7,7), S(8,8), S(9,9) ];
    end

    %% Save Space For Plotting
    g_t(:, ii) = C_t_e'*g_e;

    %extract rpy
    C_b_t = C_t_e'*C_b_e_old;
    rpy(:, ii) = extract_rpy(C_b_t);

    var_rho(:, ii) = [P(1,1), P(2,2), P(3,3)]';
    var_v_e(:, ii) = [P(4,4), P(5,5), P(6,6)]';
    var_r_e(:, ii) = [P(7,7), P(8,8), P(9,9)]';
    var_b(:, ii) = [P(10,10), P(11,11), P(12,12), P(13,13), P(14,14), P(15,15)]';
end

%% Plots
plot_idx = round(plot_until*fsamp);

figure(1);clf
subplot(3,1,1);
plot(t(1:plot_idx),(delta_th(1, 1:plot_idx)), '.');
grid on
xlabel("time (s)");
ylabel("delta Roll (rad)");
title('Delta Th for entire run.')
subplot(3,1,2);
plot(t(1:plot_idx),(delta_th(2, 1:plot_idx)), '.');
grid on
xlabel("time (s)");
ylabel("delta Pitch (rad)");
subplot(3,1,3);
plot(t(1:plot_idx), (delta_th(3, 1:plot_idx)), '.');
grid on
xlabel("t (s)");
ylabel("delta Yaw (rad)");

figure(2);clf
subplot(3,1,1);
plot(t(1:plot_idx),(delta_v(1, 1:plot_idx)), '.');
grid on
xlabel("time (s)");
ylabel("delta vel_x (m/s)");
title('Delta Vel for entire run.')
subplot(3,1,2);
plot(t(1:plot_idx),(delta_v(2, 1:plot_idx)), '.');
grid on
xlabel("time (s)");
ylabel("delta vel_y (m/s)");
subplot(3,1,3);
plot(t(1:plot_idx), (delta_v(3, 1:plot_idx)), '.');
grid on
xlabel("t (s)");
ylabel("delta vel_z (m/s)");

figure(3); clf;
subplot(3, 1, 1);
plot(t(1:plot_idx), g_t(1, 1:plot_idx), '.');
grid on;
xlabel("time (s)");
ylabel("g_x (m/s^2)");
title("Gravity in tangent frame for entire run.")
subplot(3, 1, 2);
plot(t(1:plot_idx), g_t(2, 1:plot_idx), '.');
grid on;
xlabel("time (s)");
ylabel("g_y (m/s^2)");
subplot(3, 1, 3);
plot(t(1:plot_idx), g_t(3, 1:plot_idx), '.');
grid on;
xlabel("time (s)");
ylabel("g_z (m/s^2)");

figure(4); clf;
subplot(3, 1, 1);
plot(t(1:plot_idx), r_e(1, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx),(r_e(1, 1:plot_idx)) - sqrt(var_r_e(1, 1:plot_idx)), 'k.');
plot(t(1:plot_idx),(r_e(1, 1:plot_idx)) + sqrt(var_r_e(1, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("r_x (m)");
title("Position in Earth frame for entire run.")
subplot(3, 1, 2);
plot(t(1:plot_idx), r_e(2, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx),(r_e(2, 1:plot_idx)) - sqrt(var_r_e(2, 1:plot_idx)), 'k.');
plot(t(1:plot_idx),(r_e(2, 1:plot_idx)) + sqrt(var_r_e(2, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("r_y (m)");
subplot(3, 1, 3);
plot(t(1:plot_idx), r_e(3, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx),(r_e(3, 1:plot_idx)) - sqrt(var_r_e(3, 1:plot_idx)), 'k.');
plot(t(1:plot_idx),(r_e(3, 1:plot_idx)) + sqrt(var_r_e(3, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("r_y (m)");

figure(5); clf;
subplot(3, 1, 1);
plot(t(1:plot_idx), v_e(1, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx),(v_e(1, 1:plot_idx)) - sqrt(var_v_e(1, 1:plot_idx)), 'k.');
plot(t(1:plot_idx),(v_e(1, 1:plot_idx)) + sqrt(var_v_e(1, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("v_x (m/s)");
title("Velocity in Earth frame for entire run.")
subplot(3, 1, 2);
plot(t(1:plot_idx), v_e(2, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx),(v_e(2, 1:plot_idx)) - sqrt(var_v_e(2, 1:plot_idx)), 'k.');
plot(t(1:plot_idx),(v_e(2, 1:plot_idx)) + sqrt(var_v_e(2, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("v_y (m/s)");
subplot(3, 1, 3);
plot(t(1:plot_idx), v_e(3, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx),(v_e(3, 1:plot_idx)) - sqrt(var_v_e(3, 1:plot_idx)), 'k.');
plot(t(1:plot_idx),(v_e(3, 1:plot_idx)) + sqrt(var_v_e(3, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("v_z (m/s)");

figure(6); clf;
subplot(3, 1, 1);
plot(t(1:plot_idx), rpy(1, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx),(rpy(1, 1:plot_idx)) - sqrt(var_rho(1, 1:plot_idx)), 'k.');
plot(t(1:plot_idx),(rpy(1, 1:plot_idx)) + sqrt(var_rho(1, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("roll (rad)");
title("Roll, Pitch and Yaw from body to local tangent frame for entire run.")
subplot(3, 1, 2);
plot(t(1:plot_idx), rpy(2, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx),(rpy(2, 1:plot_idx)) - sqrt(var_rho(2, 1:plot_idx)), 'k.');
plot(t(1:plot_idx),(rpy(2, 1:plot_idx)) + sqrt(var_rho(2, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("pitch (rad)");
subplot(3, 1, 3);
plot(t(1:plot_idx), rpy(3, 1:plot_idx), '.');
grid on;
hold on;
plot(t(1:plot_idx), (rpy(3, 1:plot_idx)) - sqrt(var_rho(3, 1:plot_idx)), 'k.');
plot(t(1:plot_idx), (rpy(3, 1:plot_idx)) + sqrt(var_rho(3, 1:plot_idx)), 'k.');
xlabel("time (s)");
ylabel("yaw (rad)");

figure(7);clf
subplot(3,1,1);
plot(t(1:plot_idx), b_hat(1, 1:plot_idx), '.');
grid on
hold on
plot(t(1:plot_idx), b_hat(1, 1:plot_idx) - sqrt(var_b(1, 1:plot_idx)), '.k');
plot(t(1:plot_idx), b_hat(1, 1:plot_idx) + sqrt(var_b(1, 1:plot_idx)), '.k');
title("Bias Accel");
xlabel("time (s)");
ylabel("b f_x (m/s^2)");
subplot(3,1,2);
plot(t(1:plot_idx), b_hat(2, 1:plot_idx), '.');
grid on
hold on
plot(t(1:plot_idx), b_hat(2, 1:plot_idx) - sqrt(var_b(2, 1:plot_idx)), '.k');
plot(t(1:plot_idx), b_hat(2, 1:plot_idx) + sqrt(var_b(2, 1:plot_idx)), '.k');
xlabel("time (s)");
ylabel("b f_y (m/s^2)");
subplot(3,1,3);
plot(t(1:plot_idx), b_hat(3, 1:plot_idx), '.');
grid on
hold on
plot(t(1:plot_idx), b_hat(3, 1:plot_idx) - sqrt(var_b(3, 1:plot_idx)), '.k');
plot(t(1:plot_idx), b_hat(3, 1:plot_idx) + sqrt(var_b(3, 1:plot_idx)), '.k');
xlabel("t (s)");
ylabel("b f_z (m/s^2)");

figure(8);clf
subplot(3,1,1);
plot(t(1:plot_idx), b_hat(4, 1:plot_idx), '.');
grid on
hold on
plot(t(1:plot_idx), b_hat(4, 1:plot_idx) - sqrt(var_b(4, 1:plot_idx)), '.k');
plot(t(1:plot_idx), b_hat(4, 1:plot_idx) + sqrt(var_b(4, 1:plot_idx)), '.k');
title("Bias Gyro");
xlabel("time (s)");
ylabel("b g_x (rad/s)");
subplot(3,1,2);
plot(t(1:plot_idx), b_hat(5, 1:plot_idx), '.');
grid on
hold on
plot(t(1:plot_idx), b_hat(5, 1:plot_idx) - sqrt(var_b(5, 1:plot_idx)), '.k');
plot(t(1:plot_idx), b_hat(5, 1:plot_idx) + sqrt(var_b(5, 1:plot_idx)), '.k');
xlabel("time (s)");
ylabel("b g_y (rad/s)");
subplot(3,1,3);
plot(t(1:plot_idx), b_hat(6, 1:plot_idx), '.');
grid on
hold on
plot(t(1:plot_idx), b_hat(6, 1:plot_idx) - sqrt(var_b(6, 1:plot_idx)), '.k');
plot(t(1:plot_idx), b_hat(6, 1:plot_idx) + sqrt(var_b(6, 1:plot_idx)), '.k');
xlabel("t (s)");
ylabel("b g_z (rad/s)");

figure(9); clf;
subplot(3, 1, 1);
plot(t(1:tot_samps), resid(1, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(1, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(1, 1:tot_samps)), 'k.');
title("Residual and Covariance For Attitude");
xlabel("time");
ylabel("Error Attitude X");
subplot(3, 1, 2);
plot(t(1:tot_samps), resid(2, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(2, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(2, 1:tot_samps)), 'k.');
title("Error Attitude Y");
xlabel("time");
ylabel("Error Attitude Y");
subplot(3, 1, 3);
plot(t(1:tot_samps), resid(3, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(3, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(3, 1:tot_samps)), 'k.');
title("Error Attitude Z");
xlabel("time");
ylabel("Error Attitude Z");

figure(10); clf;
subplot(3, 1, 1);
plot(t(1:tot_samps), resid(4, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(4, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(4, 1:tot_samps)), 'k.');
title("Residual and Covariance For Velocity");
xlabel("time");
ylabel("Error Velocity X [m/s]");
subplot(3, 1, 2);
plot(t(1:tot_samps), resid(5, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(5, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(5, 1:tot_samps)), 'k.');
xlabel("time");
ylabel("Error Velocity Y [m/s]");
subplot(3, 1, 3);
plot(t(1:tot_samps), resid(6, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(6, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(6, 1:tot_samps)), 'k.');
xlabel("time");
ylabel("Error Velocity Z [m/s]");

figure(11); clf;
subplot(3, 1, 1);
plot(t(1:tot_samps), resid(7, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(7, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(7, 1:tot_samps)), 'k.');
title("Residual and Covariance For Position");
xlabel("time");
ylabel("Error Position X [m]");
subplot(3, 1, 2);
plot(t(1:tot_samps), resid(8, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(8, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(8, 1:tot_samps)), 'k.');
xlabel("time");
ylabel("Error Position Y [m]");
subplot(3, 1, 3);
plot(t(1:tot_samps), resid(9, 1:tot_samps), 'b.');
grid on
hold on
plot(t(1:tot_samps), -sqrt(var_resid(9, 1:tot_samps)), 'k.');
plot(t(1:tot_samps), sqrt(var_resid(9, 1:tot_samps)), 'k.');
xlabel("time");
ylabel("Error Position Z [m]");


%% functions
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

function C_t_e = t2e(lat, lon)
    lat = lat * (pi/180);
    lon = lon * (pi/180);

    C11 = -sin(lat)*cos(lon);
    C12 = -sin(lon);
    C13 = -cos(lat)*cos(lon);
    C21 = -sin(lat)*sin(lon);
    C22 = cos(lon);
    C23 = -cos(lat)*sin(lon);
    C31 = cos(lat);
    C32 = 0;
    C33 = -sin(lat);

    C_t_e = [[C11, C12, C13];
                [C21, C22, C23];
                [C31, C32, C33]];
end

function C_b_e_new = attitude_update(C_b_e_old, alpha, tau)
    omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s) 
    wie = [0, 0, omega_ie]';

    Omega_ie_e = Skew_symmetric(wie);
    Omega_ib_b = Skew_symmetric(alpha);

    C_b_e_new = C_b_e_old*(eye(3) + Omega_ib_b) - Omega_ie_e*C_b_e_old*tau;
end

function nu_e = int_specific_force_rot(nu_b, C_b_e_old, C_b_e_new)
    nu_e = (1/2)*(C_b_e_old + C_b_e_new)*nu_b;
end

function v_e_new = velocity_update(v_e_old, nu_e, g_e, tau)
    omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s) 
    wie = [0, 0, omega_ie]';

    Omega_ie_e = Skew_symmetric(wie);

    v_e_new = v_e_old + nu_e + (g_e - 2*Omega_ie_e*v_e_old)*tau;
end

function r_e_new = position_update(r_e_old, v_e_old, v_e_new, tau)
    r_e_new = r_e_old + (v_e_old + v_e_new)*(tau/2);
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


    rpy = [(phi), (theta), (psi)]';
end

function Phi = calcStateTransition(C_b_e, alpha, r_e, g_e, lat, tau)
    Phi = eye(15);

    F_21 = -Skew_symmetric((C_b_e*alpha))*tau;

    lat = lat*(pi/180);
    e = 0.0818191908425; %WSG84 Ecentricity
    Ro = 6378137; %[m] WGS84 Equitorial Radius of Earth
    Re = Ro / sqrt(1 - e^2*sin(lat)^2);
    r_es_e = Re*sqrt(cos(lat)^2 + (1 - e^2)^2*sin(lat)^2);

    F_23 = -((2*g_e)./(r_es_e))*((r_e')./norm(r_e));

    omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s) 
    wie = [0, 0, omega_ie]';

    Omega_ie_e = Skew_symmetric(wie);

    Phi(1:3, 1:3) = eye(3) - Omega_ie_e*tau;
    Phi(1:3, 13:15) = -C_b_e*tau;
    Phi(4:6, 1:3) = F_21*tau;
    Phi(4:6, 4:6) = eye(3)-2*Omega_ie_e*tau;
    Phi(4:6, 7:9) = F_23*tau;
    Phi(4:6, 10:12) = -C_b_e*tau;
    Phi(7:9, 4:6) = eye(3)*tau;
end

function Q_d = calcINSNoise(tau, Sra, Srg, Sba, Sbg)
    Q_d = zeros(15, 15);

    Q_d(1:3, 1:3) = Srg*eye(3);
    Q_d(4:6, 4:6) = Sra*eye(3);
    Q_d(10:12, 10:12) = Sba*eye(3);
    Q_d(13:15, 13:15) = Sbg*eye(3);
    Q_d = Q_d*tau;
end

function x = inv_skew_symmetric(R)
    x = zeros(3, 1);

    x(1) = mean([-R(2, 3),R(3, 2)]);
    x(2) = mean([-R(3, 1), R(1, 3)]);
    x(3) = mean([-R(1, 2), R(2, 1)]);
end

function r_e = LLA_to_ECEF(lla) 
    lat = lla(1);
    lon = lla(2);
    alt = lla(3);

    a = 6378137;     % semi-major axis of the earth in meters (WGS 84)
    f = 1/298.257223563; % flattening factor of wgs84 (oblate) spheriod of Earth 
    e_sq = (2*f - f^2); % eccentricity squared
    Rn = a./sqrt(1-e_sq*sind(lat).^2);
    h = alt;
    x = (Rn+h).*cosd(lat).*cosd(lon);
    y = (Rn+h).*cosd(lat).*sind(lon);
    z = ((1-e_sq)*Rn + h).*sind(lat);

    r_e = [x, y, z]';
end

function lla = ECEF_to_LLA(r_e)
    x = r_e(1);
    y = r_e(2);
    z = r_e(3);

    a=6378137.0;
        % e=0.08181919092890624;
        f = 1/298.257223563; % flattening factor of wgs84 (oblate) spheriod of Earth 
        e_sq = (2*f - f^2);
        e = sqrt(e_sq);

        p=sqrt(x.^2+y.^2);
        [n,m] = size(x);
        % Initializes Newton-Raphson
        k_new=1/(1-e^2)*ones(n,m);
        err_threshold=0.0001*ones(n,m);
        err=1*ones(n,m);
        % Iterates Newton-Raphson
        while any(err>err_threshold)
            k_old=k_new;
            ci=(p.^2+(1-e^2)*z.^2.*k_old.^2).^(3/2)/(a*e^2);
            k_new=1+(p.^2+(1-e^2)*z.^2.*k_old.^3)./(ci-p.^2);
            err=abs(k_new-k_old);
        end
        k=k_new;

        lon=atan2(y,x); % Calculate longitude
        lat=atan2(z.*k,p);   % Calculate latitude
        % if lon>pi
        %     lon=lon-2*pi;   % Ensures longitude is [-pi,pi)
        % end

        Rn = a./sqrt(1-e^2*sin(lat).^2);
        alt = p./cos(lat) - Rn;

        lla = [[lat, lon]*(180/pi), alt];
end

function g_e = Gravity_ECEF(r_e)
    R_0 = 6378137; %WGS84 Equatorial radius in meters
    mu = 3.986004418E14; %WGS84 Earth gravitational constant (m^3 s^-2)
    J_2 = 1.082627E-3; %WGS84 Earth's second gravitational constant
    omega_ie = 7.292115E-5;  % Earth rotation rate (rad/s)

    % Calculate distance from center of the Earth
    mag_r = sqrt(r_e' * r_e);
    
    % If the input position is 0,0,0, produce a dummy output
    if mag_r==0
        g_e = [0;0;0];
        
    % Calculate gravitational acceleration using (2.142)
    else
        z_scale = 5 * (r_e(3) / mag_r)^2;
        gamma = -mu / mag_r^3 *(r_e + 1.5 * J_2 * (R_0 / mag_r)^2 *...
            [(1 - z_scale) * r_e(1); (1 - z_scale) * r_e(2);...
            (3 - z_scale) * r_e(3)]);
    
        % Add centripetal acceleration using (2.133)
        g_e(1:2,1) = gamma(1:2) + omega_ie^2 * r_e(1:2);
        g_e(3) = gamma(3);
        
    end
end