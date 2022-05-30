function gravity = gravitational_forces(x, y, z)
    Rp = 6356752.31425;
    Ro = 6378137.0;
    e = sqrt(1-((Rp^2)/(Ro^2)));
    wie = 7.29e-5;

    omega_ie = [[0, -wie, 0],
                [wie, 0, 0],
                [0, 0, 0]];

    [lat, lon] = lat_lon_conv(x, y, z);

    go = 9.7803253359 * ((1 + 0.001931853*sin(lat)^2)/(sqrt((1-e^2)*(sin(lat)^2))));
    gamma_o = go + omega_ie*omega_ie*[x,y,z]';
    RE = Ro / (sqrt((1-e^2)*(sin(lat)^2)));
    res = RE * sqrt(cos(lat)^2+(1-e^2)*(sin(lat)^2));
    hb = (z/sin(lat))-(1-e^2)*RE;
    gamma_ib = ((res)'/(res+hb)^2)*gamma_o;
    intermediate = [[1, 0, 0],
                    [0, 1, 0],
                    [0, 0, 0]];
    gravity = gamma_ib + omega_ie*intermediate*[x, y, z]';
end