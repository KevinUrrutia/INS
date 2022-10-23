function NED2ECEF_rot = ecefRotation(L, lam)
    R11 = -sin(L)*cos(lam);
    R12 = -sin(lam);
    R13 = -cos(L)*cos(lam);
    R21 = -sin(L)*sin(lam);
    R22 = cos(lam);
    R23 = -cos(L)*sin(lam);
    R31 = cos(L);
    R32 = 0;
    R33 = -sin(L);

    NED2ECEF_rot = [[R11, R12, R13],
                    [R21, R22, R23],
                    [R31, R32, R33]];
end