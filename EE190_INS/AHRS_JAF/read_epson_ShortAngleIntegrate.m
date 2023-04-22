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
dang= IMU(ind,6:8)*D2R;   % rad
L = length(tm);           % 

% define the stationary times; test Euler computation at those times
% initialize the rotation matrix
i_stationary = Fs:Fs:L;      % indices that IMU is stationary
T_stationary = i_stationary*T;     % times that the IMU is stationary
Ls = length(T_stationary);

%% test Euler computation
% At each stationary time, use the accel meas to compute the Euler angles
% Use teh Euler angles to compute rho. Compute the variation in each. 
I = eye(3,3);
E  = zeros(3,Ls);     % preallocate
for i=1:Ls
    E(1:3,i) = Euler_stationary_init(dv(i_stationary(i),:));   % radians
    Rn2b = calc_Rn2b(E(:,i));
    rho(:,i) = AssessAngleError(Rn2b', I);
end
Euler_deg = E*180/pi
stdE = std(E')              % std of RPY
sigE = mean(stdE(1:2))      % average RP std to get a ballpark value.
sigrho = std(rho')          % std of rho
clear rho

% use stationary interval at the start to initialize 
Na = 100;                   % number of points to average
[bg,Pg,ba,Pa,Rn2b] = InitBiasAng(dv(1:Na,:),dang(1:Na,:),grav_n*T);
Px = zeros(6,6);    
Px(1:3,1:3) = (sigE^2)*eye(3,3); % rads^2
Px(4:6,4:6) = Pg;               % (rad/sec)^2
plot_IMU(tm, dv, dang, i_stationary, L, T, bg, ba)
dang_hat = dang(:,1:3)-bg;  % bias should be near zero.
Cb2n = Rn2b';
% integrate without any corrections while stationary
[tf, xf, stdx] = AngleIntegrate(Px, dang_hat', Cb2n, T);
Lt = length(tf);
for i=1:Lt
   E = xf(:,i);
   Rn2b = calc_Rn2b(E);
   rho(:,i) = AssessAngleError(Rn2b, Cb2n');
   % TODO: Plots of rho should start at 0, but they do not. Something is
   % wrong. 
   % TODO: Use rho to correct Rn2b.  If correct, it should result in Cb2n'.
   % Ifnot, there is a bug here.
   Pc = skew(rho(:, i));
   C_corr = eye(3) + Pc;
   C_new = C_corr * Rn2b;

%    tf = isequal(C_new, Cb2n');
%    if(tf)
%        disp("Correction works correctly")
%    end
end
plot_angles(tf,xf,stdx,bg,tf,rho);




function [] = plot_angles(t, x, stdx, b, ts, rho)


figure(2)
clf
cntax =  1;
ax(cntax) = subplot(321);
plot(ts,rho(cntax,:),'.',t,stdx(cntax,:),'k',t,-stdx(cntax,:),'k')
grid on
xlabel('Time, t, sec')
ylabel('\rho_1, rad')
title('INS angles')
cntax = cntax + 1;
ax(cntax) = subplot(323);
plot(ts,rho(cntax,:),'.',t,stdx(cntax,:),'k',t,-stdx(cntax,:),'k')
grid on
xlabel('Time, t, sec')
ylabel('\rho_2, rad')
cntax = cntax + 1;
ax(cntax) = subplot(325);
plot(ts,rho(cntax,:),'.',t,stdx(cntax,:),'k',t,-stdx(cntax,:),'k')
grid on
xlabel('Time, t, sec')
ylabel('\rho_3, rad')

cntax = cntax + 1;
ax(cntax) = subplot(322);
plot([t(1),t(end)],[b(1),b(1)],t,stdx(cntax,:),'k',t,b(1)-stdx(cntax,:),'k')
grid on
xlabel('Time, t, sec')
ylabel('Bias x, rad/s')
cntax = cntax + 1;
ax(cntax) = subplot(324);
plot([t(1),t(end)],[b(2),b(2)],t,stdx(cntax,:),'k',t,b(2)-stdx(cntax,:),'k')
grid on
xlabel('Time, t, sec')
ylabel('Bias y, rad/s')
cntax = cntax + 1;
ax(cntax) = subplot(326);
plot([t(1),t(end)],[b(3),b(3)],t,stdx(cntax,:),'k',t,b(3)-stdx(cntax,:),'k')
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

% plot raw IMU data. Written as a function fo keep main clean
function [] = plot_IMU(t, dv, dang, i_stationary,L,T,bg,ba)
figure(1)
clf
cntax = 1;
ax(cntax) = subplot(321);
plot(t,dv(:,1),'.',t(i_stationary),dv(i_stationary,1),'*',t,dv(:,1)-ba(1),'.')
grid on
xlabel('Time, t, sec')
ylabel('dv_x, m/s')
title('Raw IMU - delta vel')

cntax = cntax + 1;
ax(cntax) = subplot(323);
plot(t,dv(:,2),'.',t(i_stationary),dv(i_stationary,2),'*',t,dv(:,2)-ba(2),'.')
grid on
xlabel('Time, t, sec')
ylabel('dv_y, m/s')

cntax = cntax + 1;
ax(cntax) = subplot(325);
plot(t,dv(:,3),'.',t(i_stationary),dv(i_stationary,3),'*',t,dv(:,3)-ba(3),'.')
grid on
xlabel('Time, t, sec')
ylabel('dv_z, m/s')

cntax = cntax + 1;
ax(cntax) = subplot(322);
plot(t,dang(:,1),'.',t(i_stationary),dang(i_stationary,1),'*',t,dang(:,1)-bg(1),'.')
grid on
xlabel('Time, t, sec')
ylabel('da_x, rad')
title('Raw IMU - delta angle')
legend('meas','stat','hat','Location','best')

cntax = cntax + 1;
ax(cntax) = subplot(324);
plot(t,dang(:,2),'.',t(i_stationary),dang(i_stationary,2),'*',t,dang(:,2)-bg(2),'.')
grid on
xlabel('Time, t, sec')
ylabel('da_y, rad')

cntax = cntax + 1;
ax(cntax) = subplot(326);
plot(t,dang(:,3),'.',t(i_stationary),dang(i_stationary,3),'*',t,dang(:,3)-bg(3),'.')
grid on
xlabel('Time, t, sec')
ylabel('da_z, rad')

linkaxes(ax,'x')
xlim([0,L*T])       % changes all axes to fill the horizontal space.
