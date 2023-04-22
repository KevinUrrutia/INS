function [E] = C2Euler(Cb2n)
% Aided Navigation eqns (2.45-2.47)

E = [  atan2(Cb2n(3,2), Cb2n(3,3))
       -atan(Cb2n(3,1)/sqrt(1-Cb2n(3,1)^2)) %KUA changed positions of 1,2 to what they are now
       atan2(Cb2n(2,1), Cb2n(1,1))];