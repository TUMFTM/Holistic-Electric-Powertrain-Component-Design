function dc = create_acc_cycle(t_acc, v_max, corner_speed)
% default parameters:
% t_acc = 5s for 0-100km/h
% v_max = 200 km/h
% % a2 = 1 m/s^2 acceleration from 100 to vmax minus hyperbel
    
    sample_time = 1;        % in seconds
    acc2_seconds = t_acc;       % curvature scaling (half acceleration after 100km/h)
    
    v1 = corner_speed / 3.6;
    a1 = v1 / t_acc;
    v_max = v_max / 3.6;
        
    speed_array = [0];
    v_new = 0;
    for t = 1:sample_time:1000
        if v_new < v1
            v_new = a1 * t;
        else
            acc_fact = acc2_seconds / (acc2_seconds + t) ;
            v_new = speed_array(end) + a1 * sample_time * acc_fact - 0.003 * speed_array(end);
        end

        if v_new >= v_max
            v_new = v_max;
            speed_array = [speed_array, v_new];
            break;
        end

        speed_array = [speed_array, v_new];
    end
    
    time_array = 0:sample_time:t;
    dc = {};
    dc.time = time_array';
    dc.speed = speed_array';

end