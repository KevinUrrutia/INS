function new_pos = position_update(prev_pos, prev_vel, new_vel, dt)
    new_pos = prev_pos + (prev_vel + new_vel)*(dt/2);
end