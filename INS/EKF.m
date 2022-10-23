function [x_EKF, P_mag] = EKF(delta_x, pos, vel, z, Q, P, R, dt, C, f, lat)
%     pos = x_state(begin:3);
%     vel = x_state(3:end);
    x_state = vertcat(pos, vel);

    %Calculate the Phi matrix to update the positioning
    
    %F_21 is a function of the current attitude and the current specifica
    %force estimate
    skew_elem = -(C*f);
    F_21 = [[0, -skew_elem(3), skew_elem(2)],
            [skew_elem(3), 0, -skew_elem(1)],
            [-skew_elem(2), skew_elem(1), 0]];

    %F_23 is a function of the graviational force imposed by the earth at
    %the latitude position and the current position as defined in ECEF
    %frame
    Rp = 6356752.31425;
    Ro = 6378137.0;
    w_ie = 7.29e-5; %rad/s
% 
    e = sqrt(1 - ((Rp^2)/(Ro^2)));
%     Re = Ro / sqrt(1-(e^2)*sin(lat)^2);
% 
%     hb = (pos(3) / sin(lat)) - (1-(e^2))*Re*lat;

%     g_o = 9.7803253359 * ((1 + 0.001931853*sin(lat)^2)/sqrt(1-(e^2)*sin(lat)^2));
%     g_o = [0, 0, g_o]';
    omega_ie = [[0 -w_ie, 0],
               [w_ie, 0, 0],
               [0, 0, 0]];
%     gamma_o = g_o + omega_ie*omega_ie*pos;
%     gamma = ((pos)' ./ (pos + hb).^2)*gamma_o;
% 
%     ab_pos = pos'/ abs(pos');
%     F_23 = (-(2*gamma)/pos) .* ab_pos; 
    
    geocentric_radius = Ro / sqrt(1 - (e * sin(lat))^2) *...
        sqrt(cos(lat)^2 + (1 - e^2)^2 * sin(lat)^2);
    F_23 = -2 * Gravity_ECEF(pos) /...
    geocentric_radius * (pos' / (pos' * pos));

    dt = 1;
    Phi = [[eye(3)-omega_ie*dt, zeros(3), zeros(3), zeros(3), C*dt],
          [F_21*dt, eye(3)-2*omega_ie*dt, F_23*dt, C*dt, zeros(3)],
          [zeros(3), eye(3)*dt, eye(3), zeros(3), zeros(3)],
          [zeros(3), zeros(3), zeros(3), eye(3), zeros(3)],
          [zeros(3), zeros(3), zeros(3), zeros(3), eye(3)]];

   %generate the prediction of the covariance
   Q = Phi * Q * Phi';
   P = Phi*P*Phi' + Q;

   %gernerate the measurement matrix
   H = [[zeros(3), zeros(3), -eye(3), zeros(3), zeros(3)], 
        [zeros(3), -eye(3), zeros(3), zeros(3), zeros(3)]];

   %z is the measurement of the velocity and position gotten from the gps
   %create the residual y
   y = z-x_state;

   S = H*P*H' + R; 
   K = P*H'*inv(S);
   x_EKF = delta_x + K*y;
   P_mag = (eye(15)-K*H)*P;

end