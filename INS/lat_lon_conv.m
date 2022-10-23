function [L, lam] = lat_lon_conv(x, y, z)
    Rp = 6356752.31425;
    Ro = 6378137.0;
    e = sqrt(1-((Rp^2)/(Ro^2)));

    L = atan2(z,(((1-e^2)*sqrt(x^2+y^2))));
    lam = atan2(y,x);
end