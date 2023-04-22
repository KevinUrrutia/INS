function [rho] = AssessAngleError(Rc, Rt)
% Inputs
%   - Rc is the computed rotation matrix to be assessed
%   - Rt is the "true" rotation matrix

Re = Rc' * Rt;  % = I - [rho x] 
S  = eye(3,3) - Re;
rho(1,1) = (S(3,2)-S(2,3))/2;
rho(2,1) = (S(1,3)-S(3,1))/2;
rho(3,1) = (S(2,1)-S(1,2))/2;

