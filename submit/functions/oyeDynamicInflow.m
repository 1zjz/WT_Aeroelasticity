function[vind_2, vint_2] = oyeDynamicInflow(vind, C_T1, C_T2, vint, U0, R, r,dt,glc_us)
%determine the time derivative of the induction at the annulli
    %vint is an intermediate value
    %calculate quasi-steady induction velocity
    vqst1 = -ainduction(-C_T1,0)*U0;
    %calculate current induction factor
    a = -vind/U0;
    %calculate time scales of the model
    t1 = 1.1/(1-1.3*a)*R/U0;
    t2 = (0.39-0.26*(r/R)^2)*t1;
    %calculate next-time-step quasi-steady induction velocity
    vqst2 = -ainduction(-C_T2,0)*U0;
    %calculate time derivative of intermediate velocity
    dvint_dt = (vqst1 + (vqst2-vqst1)/dt*0.6*t1 - vint)/t1;
    %calculate new intermediate velocity
    vint_2 = vint + dvint_dt*dt;
    %calculate time derivaive of the induced velocity
    dvind_dt = ((vint + vint_2)/2 - vind)/t2;
    %calculate new induced velocity
    vind_2 = vind +dvind_dt*dt;
end

    
        
    

    
    
    
    

  
    