function []=read_epson_ShortAngleIntegrate()
clear
format compact

% define constants
file = 'Epson_G370_20230407_030457.txt';
Fs =125;                        % sampling frequency, Hz
T = 1/Fs;                       % sampling period, sec
D2R = pi/180;
grav_n = [0;0;9.78];            % gravity in N-frame, m/s/s
T_stop = 22.5;                  % upper limit on time, sec

% Read data file
IMU = readmatrix(file);
L = length(IMU);
ind = 2:L;            % remove the first point
tm  = IMU(ind,1);     % sec
t = (1:L) * T;      % replace unreadable time vector
[~,im] = min(abs(t-T_stop));
ind = 1:im;
tm  = t(ind);     % sec
cnt = IMU(ind,2);
dv  = IMU(ind,3:5);   % m/s
da_nb_b= IMU(ind,6:8)*D2R;   % rad, delta ang of b wrt n in b, da_nb_b
L = length(tm);           % 

% define the stationary times; test Euler computation at those times
% initialize the rotation matrix
i_stationary = Fs:Fs:L;      % indices that IMU is stationary
T_stationary = i_stationary*T;     % times that the IMU is stationary
Ls = length(T_stationary);

% Find correct rotation: Euler computation
f = mean(dv);
E = Euler_stationary_init(f);   % radians
Rn2b_t = calc_Rn2b(E);           % true value
sigE = 8.4e-4;                   % radians
R    = sigE^2*eye(3,3);          % measurement noise covariance

% use stationary interval at the start to initialize
Na = 100; % number of points to average
[bg,Pg,ba,Pa,Rn2b] = InitBiasAng(dv(1:Na,:),da_nb_b(1:Na,:),grav_n*T);
Px = zeros(6,6);    
Px(1:3,1:3) = (sigE^2)*eye(3,3); % rads^2
Px(4:6,4:6) = Pg/(T^2);               % (rad/sec)^2
% Px(4:6, 4:6) = Pg/T; %KUA (rad/sec)
H = eye(3,6);
bh = zeros(3,Ls*Fs);        % hold estimate gyro bias
Eh = zeros(3,Ls*Fs);        % hold estimate Euler angles
sigrho = zeros(6,Ls*Fs);    % std of rotation error
for i=1:Ls
    % time propagation
    indt = (i-1)*Fs+1:i*Fs;     % one second worth of data
    dang_hat = da_nb_b(indt,:)-bg; % bg does not change in the time update
    %Cb2n = Rn2b';               % initial rotation matrix this second
    [tf, xf, stdx, Px] = AngleIntegrate(Px, dang_hat', Rn2b, T);
    bh(1:3,indt) = repmat(bg',1,Fs);
    Eh(:,indt) = xf;
    sigrho(:,indt) = stdx;
    Rn2b = calc_Rn2b(xf(:,end));% final rotation matrix
    x = zeros(6,1);            %  prior error state
    % compute measured rho using true rotation, could also use gravity to
    % compute at each time. 
    rho(:,i) = AssessAngleError(Rn2b, Rn2b_t); % this is the meas. residual

    % measurement update
    K = Px * H' * inv(H * Px * H' + R);
    Px= Px - K * (H * Px);
    x = x + K * rho(:,i);
    correction = calc_Rn2b(x(1:3));
    Rn2b = Rn2b*correction; % updated rotation matrix, see eqn. (11)
    bg   = bg + x(4:6,1)'*T; % bg is rads per T seconds, s is rad/s
    % TODO: bg is currently unstable
end

figure(5)
clf
ind  = 1:Ls*Fs;
subplot(311)
plot(tm(ind),Eh(1,:)*180/pi,'.')
ylabel('Roll, deg')
subplot(312)
plot(tm(ind),Eh(2,:)*180/pi,'.')
ylabel('Pitch, deg')
subplot(313)
plot(tm(ind),Eh(3,:)*180/pi,'.')
ylabel('Yaw, deg')
xlabel('Time, t, sec')

figure(6)
clf
subplot(321)
plot(tm(ind),sigrho(1,:),'k',tm(ind),-sigrho(1,:),'k')
subplot(322)
plot(tm(ind),bh(1,:),'.',tm(ind),bh(1,:)+sigrho(4,:),'k',tm(ind),bh(1,:)-sigrho(4,:),'k')
ylabel('Gyro bias, b_x, rad/s')
subplot(323)
plot(tm(ind),sigrho(2,:),'k',tm(ind),-sigrho(2,:),'k')
subplot(324)
plot(tm(ind),bh(2,:),'.',tm(ind),bh(2,:)+sigrho(4,:),'k',tm(ind),bh(2,:)-sigrho(4,:),'k')
ylabel('Gyro bias, b_y, rad/s')
subplot(325)
plot(tm(ind),sigrho(3,:),'k',tm(ind),-sigrho(3,:),'k')
xlabel('Time, t, sec')
subplot(326)
plot(tm(ind),bh(3,:),'.',tm(ind),bh(3,:)+sigrho(4,:),'k',tm(ind),bh(3,:)-sigrho(4,:),'k')
ylabel('Gyro bias, b_z, rad/s')
xlabel('Time, t, sec')



function [] = plot_angles_bias_once(t, x, stdx,  b)
ts = t;
figure(2)
clf
cntax =  1;
ax(cntax) = subplot(311);
plot([t(1),t(end)],[b(1),b(1)],t,b(1)+stdx(cntax,:),'k',t,b(1)-stdx(cntax,:),'k')
grid on
xlabel('Time, t, sec')
ylabel('Bias x, rad/s')
cntax = cntax + 1;
ax(cntax) = subplot(312);
plot([t(1),t(end)],[b(2),b(2)],t,b(2)+stdx(cntax,:),'k',t,b(2)-stdx(cntax,:),'k')
grid on
xlabel('Time, t, sec')
ylabel('Bias y, rad/s')
cntax = cntax + 1;
ax(cntax) = subplot(313);
plot([t(1),t(end)],[b(3),b(3)],t,b(3)+stdx(cntax,:),'k',t,b(3)-stdx(cntax,:),'k')
grid on
xlabel('Time, t, sec')
ylabel('Bias z, rad/s')

figure(3)
clf
cntax = cntax + 1;
ax(cntax) = subplot(311);
plot(t,x(1,:),'.')
grid on
xlabel('Time, t, sec')
ylabel('Roll, rad')
title('INS angles')
cntax = cntax + 1;
ax(cntax) = subplot(312);
plot(t,x(2,:),'.')
grid on
xlabel('Time, t, sec')
ylabel('Pitch, rad')
cntax = cntax + 1;
ax(cntax) = subplot(313);
plot(t,x(3,:),'.')
grid on
xlabel('Time, t, sec')
ylabel('Yaw, rad')

linkaxes(ax,'x')
xlim([0,t(end)])       % changes all axes to fill the horizontal space.
