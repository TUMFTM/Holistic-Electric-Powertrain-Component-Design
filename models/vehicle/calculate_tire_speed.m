function [torque_tire_Nm, rpm_1min] = calculate_tire_speed(time, speed, conf)
    assert(max(speed) < 200/3.6, "Speed values out of bound. Check Driving cycle");
    
    % calculate acceleration
    acc = zeros(size(speed));
    acc(2:end) = diff(speed);
    
    m_total = conf.mass_body + conf.mass_driver ...
        + conf.mass_tires;
    
    % calculate forces
    F_g = m_total * conf.gravity_acceleration;
    F_r = conf.roll_resistance * F_g;
    F_acc = conf.rot_mass_factor * m_total * acc;
    F_air = conf.cw_air * conf.frontal_area ...
        * conf.air_density / 2 * speed.^2;
    F_wheel = F_r + F_acc + F_air;
    torque_tire_Nm = F_wheel * conf.tire_radius;
    rpm_1min = speed / (2 * conf.tire_radius * pi) * 60;
    
end