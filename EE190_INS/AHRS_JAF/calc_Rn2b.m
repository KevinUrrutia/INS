function [R_n2b] = calc_Rn2b(E)
% Computes eqn. (2.43) in Aided Navigation
sr = sin(E(1));
cr = cos(E(1));
sp = sin(E(2));
cp = cos(E(2));
sy = sin(E(3));
cy = cos(E(3));
R_n2b = [ cy*cp             sy*cp       -sp
        -sy*cr+cr*sp*sr  cy*cr+sy*sp*sr  cp*sr %KUA Changed (2, 1) from -sy*cr+cy*sp*sr to current
         sy*sr+cy*sp*cr -cy*sr+sy*sp*cr  cp*cr];