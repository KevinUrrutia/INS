% Test scrip for INS
clear
% Load in truth
load("vehc_1.mat")

sig_f_imu = 10e-2;
sig_w_imu = 10e-3;
sig_gps = 1;
sig_gps_vel = 0.0001;

% Load in IMU data
iIndIMU = 1;
time_old = imu.data.time(iIndIMU:end);
dt_imu_old = imu.data.dt(1);
%fIMU = [imu.data.aX(3000:end), imu.data.aY(3000:end), imu.data.aZ(3000:end)]';
%wIMU = [imu.data.gX(3000:end), imu.data.gY(3000:end), imu.data.gZ(3000:end)]';
fIMU_old = [imu.data.aX(iIndIMU:end), imu.data.aY(iIndIMU:end), imu.data.aZ(iIndIMU:end)]';
wIMU_old = [imu.data.gX(iIndIMU:end), imu.data.gY(iIndIMU:end), imu.data.gZ(iIndIMU:end)]';

%%% Interpolate to higher data rate
dt_new = 1/100;
timeVec_old = time_old(1):dt_imu_old:time_old(end);
timeVec_new = time_old(1):dt_new:time_old(end);

% Accelerometer Interp

fIMU_x = interp1(timeVec_old, fIMU_old(1,:), timeVec_new);
fIMU_y = interp1(timeVec_old, fIMU_old(2,:), timeVec_new);
fIMU_z = interp1(timeVec_old, fIMU_old(3,:), timeVec_new);
fIMU = [fIMU_x; fIMU_y; fIMU_z];

%  v = mvnrnd(zeros(3,1), Q);
%add noise to imu data 
w_f_imu = sig_f_imu * randn(3, 18800);
fIMU = fIMU + w_f_imu;

% Gyro interp
wIMU_x = interp1(timeVec_old, wIMU_old(1,:), timeVec_new);
wIMU_y = interp1(timeVec_old, wIMU_old(2,:), timeVec_new);
wIMU_z = interp1(timeVec_old, wIMU_old(3,:), timeVec_new);
wIMU = [wIMU_x; wIMU_y; wIMU_z];
w_w_imu = sig_w_imu * randn(3, 18800);
wIMU = wIMU + w_w_imu;

time = timeVec_new;

% Load in position, velocity, attitude
iIndState = iIndIMU;
pos_e_old = truth.data.r_eb_e(iIndState:end,:)';
vel_e_old = truth.data.v_eb_e(iIndState:end,:)';
rpy_old = truth.data.E_b_n(iIndState:end,:)';

% Interp State
pos_e_x = interp1(timeVec_old, pos_e_old(1,:), timeVec_new);
pos_e_y = interp1(timeVec_old, pos_e_old(2,:), timeVec_new);
pos_e_z = interp1(timeVec_old, pos_e_old(3,:), timeVec_new);
pos_e = [pos_e_x; pos_e_y; pos_e_z];

vel_e_x = interp1(timeVec_old, vel_e_old(1,:), timeVec_new);
vel_e_y = interp1(timeVec_old, vel_e_old(2,:), timeVec_new);
vel_e_z = interp1(timeVec_old, vel_e_old(3,:), timeVec_new);
vel_e = [vel_e_x; vel_e_y; vel_e_z];

rpy_x = interp1(timeVec_old, rpy_old(1,:), timeVec_new);
rpy_y = interp1(timeVec_old, rpy_old(2,:), timeVec_new);
rpy_z = interp1(timeVec_old, rpy_old(3,:), timeVec_new);
rpy = [rpy_x; rpy_y; rpy_z];

% Initialize
nav_frame = 'ECEF';
accel_noise_PSD = 0;
gyro_noise_PSD = 0;
ini_pos_b = pos_e(:,1);
ini_vel_b = vel_e(:,1);
ini_rpy = rpy(:,1);

%Create GPS data at 1 Hz
GPS_pos = pos_e(:, 1:100:end);
GPS_vel = vel_e(:, 1:100:end);

w_gps = sig_gps * randn(3, 188);
w_gps_vel = sig_gps_vel * randn(3, 188);
GPS_vel = GPS_vel + w_gps_vel;
GPS_pos = GPS_pos + w_gps;

%initialize the covariance matrix
% Q_attitude = blkdiag(0.00057^2, 0.00057^2, 0.00057^2);
% Q_specific_force = blkdiag(0.00011^2, 0.00011^2, 0.00011^2);
% Q_attitude = blkdiag(0, 0, 0);
% Q_specific_force = blkdiag(0, 0, 0);

% Initialize memory for saving position, velocity
numIts = 18000;
saved_pos_b = zeros(3,numIts);
saved_vel_b = zeros(3,numIts);
saved_pos_b(:,1) = ini_pos_b;
saved_vel_b(:,1) = ini_pos_b;

ini_rpy_rot = eurlerRotation(ini_rpy(3), ini_rpy(2), ini_rpy(1));
[ini_lat, ini_lon] = lat_lon_conv(ini_pos_b(1), ini_pos_b(2), ini_pos_b(3));
NED2ECEF_rot = ecefRotation(ini_lat, ini_lon);
ini_attitude = NED2ECEF_rot*ini_rpy_rot;

prev_attitude = ini_attitude;
prev_pos_ecef = ini_pos_b;
prev_vel = ini_vel_b;

%initialize the covariance matrices 
% Q = [[0.00057*eye(3), zeros(3), zeros(3), zeros(3), zeros(3)], 
%      [zeros(3), 0.00057*eye(3), zeros(3), zeros(3), zeros(3)], 
%      [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)], 
%      [zeros(3), zeros(3), zeros(3), 0.00011*eye(3), zeros(3)], 
%      [zeros(3), zeros(3), zeros(3), zeros(3), 0.00011*eye(3)]];
% S_rg = sig_w_imu^2 * dt_new;
% S_ra = sig_f_imu^2 * dt_new;
% Q = [[S_rg*eye(3), zeros(3), zeros(3), zeros(3), zeros(3)], 
%      [zeros(3), S_ra*eye(3), zeros(3), zeros(3), zeros(3)], 
%      [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)],
%      [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)], 
%      [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)]];


%initialize R as the coariance of the noise
% R = blkdiag(eye(3)*0.00057^2, eye(3)*0.00011^2);
R = blkdiag(eye(3)*sig_gps^2, eye(3)*sig_gps_vel^2);

%intialized delta state vector
prev_delta_x = zeros(15,1);

%initilaize P as the block diagnal of the variances for both the imu and
%the gnss
P = blkdiag(eye(3)*0.01, eye(3)*0.001, eye(3), eye(3)*0, eye(3)*0);

% for ii = 2:numIts
cnt = 2;
for ii = 2:numIts
    dt = time(ii) - time(ii-1);

    if(mod(ii-1, 100) ~=  0)
%       if(true)
        %take the new rpy and transform to the ECEF frame
        new_attitude = attitude_update(prev_attitude, wIMU(3, ii), wIMU(2, ii), wIMU(1, ii), dt);

        %update the specific force
        specific_force_ecef = specific_force_update(prev_attitude, new_attitude, fIMU(:, ii));

        %update the gravity
        gravity = Gravity_ECEF(prev_pos_ecef);

        %check gravity
        gravity_body = NED2ECEF_rot *ini_rpy_rot* gravity;

        %update the velocity
        new_vel = velocity_update(prev_vel, specific_force_ecef, gravity, dt);

        %update the position
        new_pos = position_update(prev_pos_ecef, prev_vel, new_vel, dt);


        % Save propagated position and velocity
        saved_pos_b(:,ii) = new_pos;
        saved_vel_b(:,ii) = new_vel;

        prev_attitude = new_attitude;
        prev_pos_ecef = new_pos;
        prev_vel = new_vel;  

    else
        %measurement model is the gnss vel and position subtracted by the
        %INS solution
        S_rg = sig_w_imu^2 * dt_new;
        S_ra = sig_f_imu^2 * dt_new;
        Q = [[S_rg*eye(3), zeros(3), zeros(3), zeros(3), zeros(3)], 
            [zeros(3), S_ra*eye(3), zeros(3), zeros(3), zeros(3)], 
            [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)],
            [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)], 
            [zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)]];
%         Q = Q * dt;
        z = vertcat(GPS_pos(:, cnt), GPS_vel(:, cnt));
        C = prev_attitude;
        f = fIMU(:, ii);
        [lat, lon] = lat_lon_conv(prev_pos_ecef(1), prev_pos_ecef(2), prev_pos_ecef(3));
        [delta_x, new_P] = EKF(prev_delta_x, prev_pos_ecef, prev_vel, z, Q, P, R, dt, C, f, lat);
        P = new_P;

        new_pos = prev_pos_ecef - delta_x(7:9);
        new_vel = prev_vel - delta_x(4:6);
        delta_euler = delta_x(1:3);
        skew_euler = [[0, -delta_euler(3), delta_euler(2)],
                      [delta_euler(3), 0, -delta_euler(1)],
                      [-delta_euler(2), delta_euler(1), 0]];

        new_attitude = (eye(3) - skew_euler)*prev_attitude;
%         prev_delta_x = delta_x;
        saved_pos_b(:,ii) = new_pos;
        saved_vel_b(:,ii) = new_vel;

        prev_attitude = new_attitude;
        prev_pos_ecef = new_pos;
        prev_vel = new_vel; 
        cnt = cnt+1;
    end
end

% Calculate errors
error_3D = sqrt(sum((saved_pos_b - pos_e(:,1:numIts)).^2, 1));
figTime = (0:numIts-1)*dt;
figure;plot(figTime, error_3D)

%Plot truth and estimate
figure;plot3(pos_e(1,1:numIts)-pos_e(1,1),pos_e(2,1:numIts)-pos_e(2,1),pos_e(3,1:numIts)-pos_e(3,1));
hold on
plot3(saved_pos_b(1,:)-saved_pos_b(1,1),saved_pos_b(2,:)-saved_pos_b(2,1),saved_pos_b(3,:)-saved_pos_b(3,1));

%plot truth and estimate in local tangent plane
%get the longitude and latitude of everypoint and transform the frame of
%each point
saved_pos_b_local = zeros(3,numIts);
pos_e_local = zeros(3,numIts);

for ii = 1:numIts
    [lat, lon] = lat_lon_conv(saved_pos_b(1,ii), saved_pos_b(2,ii), saved_pos_b(3,ii));
    NED2ECEF_rot = ecefRotation(lat, lon);
    saved_pos_b_local(:, ii) = NED2ECEF_rot' * (saved_pos_b(:, ii) - saved_pos_b(:,1));

    [lat, lon] = lat_lon_conv(pos_e(1,ii), pos_e(2,ii), pos_e(3,ii));
    NED2ECEF_rot = ecefRotation(lat, lon);
    pos_e_local(:, ii) = NED2ECEF_rot' * (pos_e(:, ii) - pos_e(:,1));


end

%plot the truth and estimate in the local tangent frame
% z = zeros(numIts);
% figure;plot3(saved_pos_b_local(1,:),saved_pos_b_local(2,:),saved_pos_b_local(3,:));
% hold on
% plot3(pos_e_local(1,:),pos_e_local(2,:),pos_e_local(3,:));

%plot the error in the local tangent frame
% error_3D_NED = sqrt(sum((saved_pos_b_local - pos_e_local(:,1:numIts)).^2, 1));
% figTime = (0:numIts-1)*dt;
% figure;plot(figTime, error_3D_NED)

% figure;plot(figTime, error_3D);
% hold on
% plot(figTime, error_3D_NED)


