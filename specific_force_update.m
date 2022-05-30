function new_specific_force = specific_force_update(prev_attitude, new_attitude, specific_force_body, Q)
    %generate the gaussian noise for the specific force update
    v = mvnrnd(zeros(3,1), Q);

    new_specific_force = (1/2) * (prev_attitude + new_attitude)*specific_force_body;
    new_specific_force = new_specific_force + v';
    
end