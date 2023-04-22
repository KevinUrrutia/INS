function [E] = Euler_stationary_init(v)
% Assuming that the IMU is stationary, this function uses the specific
% force vector to initialize the Euler angles
if v(3)<0
    v = -1*v;
end
Nv= norm(v(2:3));
E = [ atan2(v(2),v(3))
      atan(-v(1)/Nv)
      0];
