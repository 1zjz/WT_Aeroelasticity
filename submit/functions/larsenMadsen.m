function[vind_2] = larsenMadsen(vind, C_T2, U0, R,dt,glc_us)
%determines the time derivative of the induction at the annulli

    %calculate velocity wake 
    V_wake = U0 + vind;
    %calculate time scales of the model
    t1 = 0.5*R/V_wake;
    %calculate next-time-step quasi-steady induction velocity
    vqst2 = -ainduction(-C_T2,0)*U0;
    %calculate new induced velocity
    vind_2 = vind*exp(-dt/t1) + vqst2*(1-exp(-dt/t1));
end