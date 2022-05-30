function new_velocity = velocity_update(prev_vel, specific_force, gravity, dt)
%     gravity = gravitational_forces(x, y, z);
    wie = 7.29e-5;

    omega_ie = [[0 -wie 0],
                [wie 0 0],
                [0 0 0]];

%     new_velocity = prev_vel + (specific_force + gravity*([x,y,z]' - 2*omega_ie*prev_vel))*dt;
    new_velocity = prev_vel + (specific_force + gravity - 2.*omega_ie*prev_vel)*dt;
end