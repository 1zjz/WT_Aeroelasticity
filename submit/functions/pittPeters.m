function[vind_2,dvind_dt] = pittPeters(C_T, vind, U0, R, dt, glc_us)
%calculate the time derivative of the induction da/dt
    a = -vind/U0; %determine the induction coefficient for the time step {i-1}
    C_Tn = -getCT(a, 0, 0); %(a, glc, gamma) calculate the thrust coefficient 
    % from the induction for the time step {i-1}
    dvind_dt =  (C_T-C_Tn)/(16/(3*pi))*(U0^2/R); 
    %calculate the time derivative of the induction velocity
    vind_2 = vind + dvind_dt*dt %calculate the induction at time {i} by time integration

end
