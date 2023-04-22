function []=read_epson()
clear
format compact

% define constants
file = 'Epson_G370_20230407_030457.txt';
Fs =125;                        % sampling frequency, Hz
T = 1/Fs;                       % sampling period, sec
D2R = pi/180;
grav_n = [0;0;9.78];            % gravity in N-frame, m/s/s

% Read data file
IMU = readmatrix(file);
L = length(IMU);
ind = 2:L;
tm  = IMU(ind,1);     % UTC
cnt = IMU(ind,2);
dv  = IMU(ind,3:5);   % m/s
dang= IMU(ind,6:8)*D2R;   % rad
L = L-1;            % first point is discarded
t = (1:L) * T;      % replace unreadable time vector

% define the stationary times; test Euler computation at those times
% initialize the rotation matrix
i_stationary = 15*Fs:20*Fs:L;      % indices that IMU is stationary
T_stationary = i_stationary*T;  % times that the IMU is stationary
Ls = length(T_stationary);
plot_IMU(t, dv, dang, i_stationary,L,T)

% test Euler computation
E  = zeros(3,Ls);     % preallocate
for i=1:Ls
    E(1:3,i) = Euler_stationary_init(dv(i_stationary(i),:));   % radians
    Rn2b = calc_Rn2b(E(:,i));
end
Rn2b = calc_Rn2b(E(:,1));   % initialize rotation
Euler_deg = E*180/pi

% use stationary interval at the start to initialize
ind = (1:15)*Fs;    
[bg,Pg,ba,Pa,Rn2b] = InitBiasAng(dv(ind,:),dang(ind,:),grav_n*T);





% plot raw IMU data. Written as a function fo keep main clean
function [] = plot_IMU(t, dv, dang, i_stationary,L,T)
figure(1)
clf
cntax = 1;
ax(cntax) = subplot(321);
plot(t,dv(:,1),'.',t(i_stationary),dv(i_stationary,1),'*')
grid on
xlabel('Time, t, sec')
ylabel('dv_x, m/s')
title('Raw IMU - delta vel')

cntax = cntax + 1;
ax(cntax) = subplot(323);
plot(t,dv(:,2),'.',t(i_stationary),dv(i_stationary,2),'*')
grid on
xlabel('Time, t, sec')
ylabel('dv_y, m/s')

cntax = cntax + 1;
ax(cntax) = subplot(325);
plot(t,dv(:,3),'.',t(i_stationary),dv(i_stationary,3),'*')
grid on
xlabel('Time, t, sec')
ylabel('dv_z, m/s')

cntax = cntax + 1;
ax(cntax) = subplot(322);
plot(t,dang(:,1),'.',t(i_stationary),dang(i_stationary,1),'*')
grid on
xlabel('Time, t, sec')
ylabel('da_x, rad')
title('Raw IMU - delta angle')

cntax = cntax + 1;
ax(cntax) = subplot(324);
plot(t,dang(:,2),'.',t(i_stationary),dang(i_stationary,2),'*')
grid on
xlabel('Time, t, sec')
ylabel('da_y, rad')

cntax = cntax + 1;
ax(cntax) = subplot(326);
plot(t,dang(:,3),'.',t(i_stationary),dang(i_stationary,3),'*')
grid on
xlabel('Time, t, sec')
ylabel('da_z, rad')

linkaxes(ax,'x')
xlim([0,L*T])       % changes all axes to fill the horizontal space.
