function dwdt = omegaequations(t, w)
    Ix = 0.00004942;
    Iy = 0.00004317;
    Iz = 0.00000901;
    wx = w(1);
    wy = w(2);
    wz = w(3);
    phi = w(4);
    theta = w(5);
    psi = w(6);
    dwdt = [(Iy-Iz)/Ix*wy*wz; (Iz-Ix)/Iy*wz*wx; (Ix-Iy)/Iz*wy*wx;
        1/sin(theta)*(wx*sin(psi)+wy*cos(psi)); wx*cos(psi)-wy*sin(psi);
        wz-cos(theta)/sin(theta)*(wx*sin(psi)+wy*cos(psi))];