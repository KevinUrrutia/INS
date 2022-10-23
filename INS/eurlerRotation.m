function euler_rot_mat = eurlerRotation(psi, phi, theta)
    R11 = cos(psi)*cos(theta);
    R12 = cos(psi)*sin(phi)*sin(theta)-sin(psi)*cos(theta);
    R13 = cos(phi)*sin(theta)*cos(psi)+sin(psi)*sin(phi);
    R21 = sin(psi)*cos(theta);
    R22 = sin(psi)*sin(phi)*sin(theta)+cos(phi)*cos(psi);
    R23 = sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi);
    R31 = -sin(theta);
    R32 = cos(theta)*sin(phi);
    R33 = cos(phi)*cos(theta);

    euler_rot_mat = [[R11, R12, R13],
                     [R21, R22, R23],
                     [R31, R32, R33]];

end