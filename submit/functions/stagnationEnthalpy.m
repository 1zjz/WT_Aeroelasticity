function[h1, h2, h3, h4] = stagnationEnthalpy(V_inf, rho, p_inf, a, CT, A)
    h1 = p_inf/rho + V_inf^2/2;
    h2 = h1;
    
    %change in stagnational anthalpy due to external force at the actuator
    %disk: this is the thrust
        
    T = 0.2*rho*V_inf^2*A*CT;
    p3 = -T/A + p_inf;

    h3 = p3/rho + ((V_inf*(1-a))^2)/2;
    h4 = h3;

end