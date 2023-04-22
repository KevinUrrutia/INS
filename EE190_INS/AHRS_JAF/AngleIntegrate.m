function [tf, xf, stdx, Pf] = AngleIntegrate(P0, da_nb_b, Cn2b, T)
% Inputs:
%   P0 initial error state covariance
%   da change in angle, 3 x N matrix
%   Cb2n initial rotation matrix
%   T time step
%
% Outputs:
%   xf final state
%   Pf final error state covariance

[~,N] = size(da_nb_b);
stdx= zeros(6,N);
P = P0;
stdx(:,1) = sqrt(diag(P0));
Elr = zeros(3,N);
Cb2n = Cn2b';
Elr(:,1) = C2Euler(Cb2n);
for i=2:N
    [phi,Qd] = calc_ang_phi_Qd(Cb2n,T);
    P = phi * P * phi' + Qd;
    stdx(:,i) = sqrt(diag(P));
    Cn2b = ang_step(Cn2b,da_nb_b(:,i-1));    % integrate the rotation matrix
    Cb2n = Cn2b';
    Elr(:,i) = C2Euler(Cb2n);
end
xf = Elr;   % only return the Euler angles to save memory
Pf = P;
tf = (T:T:N*T);

function [phi,Qd] = calc_ang_phi_Qd(Cb2n,T)
I=eye(3,3);
Z     = zeros(3,3);
F     = [Z  -Cb2n
         Z  Z];
G     = [-Cb2n  Z
         Z     I];
Q     = [(1e-6)*I Z
         Z  1e-10*I ]  ;   % TODO: Fix. These values are made up.
Z     = 0*Q;
% Eqn. (4.113) in Aided Navigation
Chi = [-F   G*Q*G'
        Z   F']*T;
Ups = expm(Chi);
phi = Ups(7:12,7:12)';
Qd  = phi * Ups(1:6,7:12);

function [Co] = ang_step(Ci,da)
% Eqn. (2.70) in Aided Navigation for Cn2b
%   books a-frame is b-frame here
%   books b-frame is n-frame here
% da is in radians, da_nb_b
I   = eye(3,3);
Ups = v2skew(da);   %  Eqn. (2.64)
na  = norm(da); 
sa  = sin(na);
ca  = cos(na);
if na == 0
    saoa = 0;
    caoa = 0.5;
else
    saoa = sa/na;
    caoa = (1-ca)/na/na;
end
Co = (I - saoa*Ups + caoa*Ups*Ups) * Ci;

function [S] = v2skew(v)
S= [ 0 -v(3) v(1)
    v(3) 0  -v(2)
   -v(1) v(2) 0];