function new_attitude = attitude_update(prev_attitude, psi, phi, theta, dt)
   %generate the gaussian noise for the attitude update
%     v = mvnrnd(zeros(3,3), Q);

    wie = 7.29e-5;
    omega_ib = [[0 -psi phi],
                [psi, 0, -theta],
                [-phi, theta, 0]];

    omega_ie = [[0 -wie 0],
                [wie 0 0],
                [0 0 0]];

    new_attitude = prev_attitude*(eye(3) + omega_ib*dt)-omega_ie*prev_attitude*dt;
%     new_attitude = new_attitude + v;


end