function [bias_g,P_g,bias_a,P_a,Rn2b]=InitBiasAng(f,w,g)
% Inputs: 
%   f is an N x 3 matrix of specific forces or delta velocities
%   w is an N x 3 matrix of angular rates or delta angles
%   g is the 1 x 3 assumed gravity vector in navigation frame

N = length(f);
if N~= length(w)
    error("InitBiasAng: f and w must have the same dimensions.")
end
% processs gyros
[bias_g,P_g,~] = CompCheckStats(w,N,0,g);
% process accelerometers
[bias_a,P_a,Rn2b] = CompCheckStats(f,N,1,g);



function [avg,cov,Rn2b] = CompCheckStats(x,N,Type,g)
% Compute a typical value and the std of x
% When x is the specific force, also init the rotation matrix.
avg_x = mean(x);
if Type == 1    % accelerometer
    Elr = Euler_stationary_init(avg_x);
    Rn2b = calc_Rn2b(Elr);
    x = x + (Rn2b * g)';    % remove gravity
    avg_x = median(x);      % recompute median
else
    Rn2b = [];
end
res_x = x - avg_x;
std_x = std(res_x);
wrst = max(abs(res_x));
if max(wrst) > 5*std_x
    wrst, std_x
    display("InitBiasAng: large residual")
end
avg = avg_x;  % 
cov = diag(std_x.*std_x)/(N);
