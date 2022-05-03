function [cP, cT,MB_normal,MB_tangential, p_normal, p_tangential] = ...
    BIG_BEMalg(R,B,rho,V0,lambda,theta,Cl,Cd,aoa_ser,glc, r, c, beta)
    %% Description
    % This function implements the BEM algorithm
    %% General Information
    % Version: 1
    % Date: 20/09/2021
    % Authors: Sowmya, Philipp, Carlos
    % Denmark Technical University (DTU)
    % Wind Turbine Technologies and Aerodynamics
    % Assignment 1
    %% Version 2:
    % Date: 08/03/2022
    % Description: The code has been adapted for the second assignment of
    % Rotor / Wake Aerodynamics course of TU Delft
    %% Function dictionary
    % _____________________________________________________________________
    % INPUTS
    % - r       ---> Radial blade element position [m]
    % - R       ---> Radius of the rotor [m]
    % - B       ---> Number of blades [-]
    % - rho     ---> Air density [kg/m3]
    % - V       ---> Mean wind speed [m/s]
    % - omega   ---> Rotational speed [rad/s]
    % - theta   ---> Pitch angle [deg]
    % - beta    ---> Local twist angle [deg] 
    % - c       ---> Chord of the blade element airfoil [m]
    % - CL      ---> Array with lift coefficient [-]
    % - CD      ---> Array with drag coefficient [-] 
    % - CM      ---> Array with moment coefficient [-]
    % - aoa_ser ---> Series of aoa for which CL,CD, and CM are given [deg]
    % - glc     ---> Glauert correction flag [1 · else]
    % - filebem ---> 
    % _____________________________________________________________________
    % OUTPUTS
    % - a                  ---> Induced factor normal [-]
    % - aprime             ---> Induced factor tangential [-]
    % - pt                 ---> Local tangential force [N/m]
    % - pn                 ---> Local normal force [N/m]
    % - F                  ---> Correction factor [-]
    % _____________________________________________________________________
    % AUXILIARY
    % - A             ---> Auxiliary matrix for storing the data [-]
    % - count         ---> Auxiliary variable to count BEM iterations [-]
    % - eps           ---> Auxiliary variable for convergence [-]
    % - eps_lim       ---> Tolerance limit for BEM convergence [-]
    % - phi           ---> Flowangle [deg]
    % - r_e           ---> Radial elemental position [m]
    % - c_e           ---> Chord of blade element    [m]
    % - beta_e        ---> Twist of blade element    [deg]
    
    % _____________________________________________________________________
    % *********************************************************************
    % Algorithm ===========================================================
    % Step 1 - Initialize a and aprime
    % Step 2 - Compute Flow Angle $Phi$ from Equ. 6.7
    %          Extra: Compute Factor F (Prandtl)
    % Step 3 - Compute local AoA from Equ. 6.6
    % Step 4 - Read C_l($alpha$) and C_d($alpha$) from table look-up
    % Step 5 - Compute C_n and C_t from Equ. 6.12 and 6.13
    % Step 6 - Calculate a and aprime from Equ. 6.23 and 6.24
    %          Prandtl's Tip Loss Factor: Use Equ. 6.35 and 6.36
    % Step 7 - IF a and aprime has changed more than tolerance, go to Step 2,
    %          ELSE Finish
    % Step 8 - Integrate over Blades for resulting Thrust and PowerThis Algo
    % =====================================================================
    %% Parameters
    % *********************************************************************
    eps_lim = 10e-10; % Set a tolerance limit for the convergence
    Nitmax = 1000;    % Maximum number of iterations
    % *********************************************************************
    %% Operation
    % *********************************************************************
    
    % Number of blade elements
    NE = length(r);
    
    PNplot = zeros((NE+1),1);
    PTplot = zeros((NE+1),1);
    a_plot = zeros(NE,1);
    aprime_plot = zeros(NE,1);
    F_plot = zeros(NE,1);
    
    
    mbending_normal = zeros(NE+1,1);
    mbending_tangential = zeros(NE+1,1);


    for ii = 1:NE % Loop over all the blade elements
        
        % Take the ii blade position values
        r_e     = r(ii);   
        c_e     = c(ii);    
        beta_e  = beta(ii); 
        
        %% Here it starts the original BEM
        % In the memory of the good old times 
        % You are entering Legoland, enjoy

        % Calculate the rotor speed for the given tip-speed ratio
        omega = lambda*V0/R; 
        
        % Initialize variables
        a       = 0;                   % Step 1 · Induced Factor normal
        aprime  = 0;                   % Induced Factor tangential
        eps     = 1;                   % Error(must be larger than eps_lim)
        count   = 0;                   % Iteration counter
        % Step 2
        % Calculate the inflow angle
        phi     = atand(((1-a)*V0) / ((1+aprime)*omega*r_e));
        
        while eps > eps_lim && count <= Nitmax % While the convergence is 
                                               % not achieved keep going
            
%             % Step 2
%             % Calculate the inflow angle
%             phi     = atand(((1-a)*V0) / ((1+aprime)*omega*r_e));
            
            % Prandtl's Tip Loss
            % empirical value, NOT acosd. Also include abs(phi)

%             aux_t   = -B/2*(1-mu)/mu*sqrt((1+(lambda*mu)^2)/(1-a)^2);
%             F_t_test = 2/pi*acos(exp(aux_t));
%             aux_r   = -(B/2)*(mu - muR)/mu*sqrt(1 + ((lambda*mu)^2)/((1 - a)^2));
%             F_r_test     = 2/pi*acos(exp(aux_r));
%             F = F_r_test * F_t_test;
            F_t = (2/pi)*acos(exp(-(B/2)*((R-r_e)/(r_e*abs(sind(phi))))));
            F_r = (2/pi)*acos(exp(-(B/2)*((r_e - r(1))/(r_e*abs(sind(phi))))));
            F = F_r*F_t;

            

            if ~isreal(F)
                disp('It became complex')
            elseif F < 10e-3
                F = 0.001;
            end
            
            phi = atan(R/lambda*r_e*(1-a*F)/(1+aprime/F));

            % Step 3
            % Compute the local angle of attack in degrees
            aoa     = phi - beta_e - theta;
            
            % Step 4 given Cl and Cd
            sigma   = (c_e*B)/(2*pi*r_e);
            Cl_af = interp1(aoa_ser, Cl, aoa);
            Cd_af = interp1(aoa_ser, Cd, aoa);
%             Cm_af = interp1(aoa_ser, Cm, aoa);
            
            % Step 5
            Cn      = Cl_af*cosd(phi) + Cd_af*sind(phi);  %Normal force                    
            Ct      = Cl_af*sind(phi) - Cd_af*cosd(phi);  %Tangential force                  
            
            if glc == 0 || abs(count) <= eps_lim
                
                CT  = ((1-a).^2 * Cn*sigma)/((sind(phi)).^2);
                % Step 6
                an      = 1 / (( (4*(sind(phi))^2) / (sigma*Cn) ) +1);
                
            else
                
                if a <= 1/3
                    an = 1/(((4*F*(sind(phi)).^2)/(sigma*Cn))+1);
                    CT = 4*an*(1-an)*F;
                else
                    astar = CT/(4*F*(1-0.25*(5-3*a)*a));
                    an = 0.1*astar + (1-0.1)*a;
                    CT = (1-a)^2*Cn*sigma / (sind(phi))^2;
                    
                end
                
            end
            
            aprime2 = 1 / ( ( (4*F*sind(phi)*cosd(phi)) / (sigma*Ct) ) -1);
            
            eps     = max(abs(a - an),abs(aprime - aprime2));           % Step 7
            a       = an;
            aprime  = aprime2;
            
            
            count   = count + 1;
        end
        
        Vrel    = V0*(1-a) / sind(phi);                             % Step 8
        pN      = 0.5 * rho * Vrel^2 * c_e * Cn;
        pT      = 0.5 * rho * Vrel^2 * c_e * Ct;
        
        %% Here it finishes the original BEM Code
        
        PNplot(ii) = pN;
        PTplot(ii) = pT;
        a_plot(ii) = a;
        aprime_plot(ii) = aprime;
        F_plot(ii) = F;

        % Calculate the bending moments with respect to the root of the 
        % blade
        mbending_normal(ii) = PNplot(ii) * (r(ii) - r(1));
        mbending_tangential(ii) = PTplot(ii) * (r(ii) - r(1));
    end

    figure
    plot(r, F_plot)
    
    PNplot(NE+1) = 0;
    PTplot(NE+1) = 0;

    p_normal     = PNplot(1:NE);
    p_tangential = PTplot(1:NE);
    
    r_full = [r R];
    power  = lambda*V0/R*B*trapz(r_full, r_full.*PTplot');
    cP     = power/(0.5*rho*(V0^3)*pi*(R^2));
    thrust = B * trapz(r_full,PNplot);
    cT     = thrust/ (0.5*rho*V0^2*pi*R^2);
    
    % Integrate the distribution of bending moments to obtain the bending
    % moment at the root of the blade
    MB_normal     = trapz(r_full,mbending_normal);
    MB_tangential = trapz(r_full,mbending_tangential);

end