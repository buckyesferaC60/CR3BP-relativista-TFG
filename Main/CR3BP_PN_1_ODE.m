function dydt = CR3BP_PN_1_ODE(t,p,F,mu,kappa)
x = p(1);
y = p(2);
z = p(3);
vel_x = p(4);
vel_y = p(5);
vel_z = p(6);

%Campo de vectores
v_x = vel_x;%round(vel_x,15,"decimals");
v_y = vel_y;%round(vel_y,15,"decimals");
v_z = vel_z;%round(vel_z,15,"decimals");

if kappa == 0
    v_velx = F(1).func(x,y,vel_y,z);%round(F(1).func,15,"decimals");
    v_vely = F(2).func(x,vel_x,y,z);%round(F(2).func,15,"decimals");
    v_velz = F(3).func(x,y,z);%round(F(3).func,15,"decimals");
else
    v_velx = F(1).func(x,vel_x,y,vel_y,z,vel_z);%round(F(1).func,15,"decimals");
    v_vely = F(2).func(x,vel_x,y,vel_y,z,vel_z);%round(F(2).func,15,"decimals");
    v_velz = F(3).func(x,vel_x,y,vel_y,z,vel_z);%round(F(3).func,15,"decimals");
end





%fprintf(".\n") % For debugging purposes

%Output
dydt = [v_x; v_y; v_z; v_velx; v_vely; v_velz];

    %if (((x+mu)^2 + y^2 + z^2)^0.5 < 0.01 )||( ((x-1+mu)^2 + y^2 + z^2)^0.5 < 0.01)
        %dydt = 0*dydt;
    %end

end

