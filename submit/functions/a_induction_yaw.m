function [a_i_yaw] = a_induction_yaw(CT, F)
    a_c = 0.2;
    a_i_yaw = (1 - sqrt(1 - CT*F))/2;
    if a_i_yaw > a_c
        a_i_yaw = (CT - 4*a_c^2*F)/(4*F - 8*a_c*F);
    end
end
