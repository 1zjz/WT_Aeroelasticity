function [a , aline, fnorm , ftan, gamma, inflowangle, aoa, Cq] = ...
         solveStreamtube(Uinf, r1_R, r2_R, rootradius_R, tipradius_R, Omega,... 
                   Radius, NBlades, chord, twist, polar_alpha, polar_cl,...
                   polar_cd, yaw_angle, psi, pitch, UsePrandtl)
    
%     Solve balance of momentum between blade element load and loading in 
%     the streamtube
%     Input variables:
%     Uinf - wind speed at infinity
%     r1_R,r2_R - edges of blade element, in fraction of Radius ;
%     rootradius_R, tipradius_R - location of blade root and tip, in fraction of Radius ;
%     Radius is the rotor radius
%     Omega - rotational velocity
%     NBlades - number of blades in rotor
%     polar_alpha, polar_cl, polar_cd - Aerofoil characteristics
%     yaw_angle - yaw angle of the wind turbine [deg]
%     psi - azimuthal position
%     pitch - pitch angle of the blade

    % Area of the portion
    Area = pi*((r2_R*Radius)^2-(r1_R*Radius)^2); 
    r_R = (r1_R+r2_R)/2; % Centroide
    % Initialize variables
    aline = 0.0;    % Tangential induction factor
    a = 0.0;        % Axial induction factor
    N = 100;
    Erroriterations = 0.00001; % Error limit for iteration process
                               % in absolute value of induction
    yaw = yaw_angle*pi/180; % Use yaw angle in radians
    psi = psi*pi/180;       % Use azimuthal angle in radians
    
    for ii = 1:N    % Loop until convergence is achieved
        
        % Obtain the velocities at the rotor
        % Calculate skew angle based on the azimuthal position, yaw angle
        % and axial induction factor, based on Coleman's expression
        chi = (0.6*a+1)*yaw;
        K  = 2*tan(chi/2);
        % Induced non-uniform axial velocity 
        u_yaw = K*r_R*sin(psi);
        % Axial velocity at rotor
        vnorm = Uinf*(cos(yaw)-a*(1+u_yaw));
        % Tangential velocity at rotor
        vtan = (1+aline)*Omega*r_R*Radius... % Rotational component
                    + Uinf*sin(yaw)*cos(psi); % Yaw downstream component
        % Calculate the inflow angle
        inflowangle = atand(vnorm/vtan);
        % Calculate loads in blade segment in 2D (N/m)
        [fnorm, ftan, gamma, aoa] = loadBladeElement(vnorm, vtan,chord,...
                                twist, polar_alpha, polar_cl, polar_cd, pitch);

        % 3D force in axial direction
        load3Daxial =fnorm*Radius*(r2_R-r1_R)*NBlades; 

        % Calculate new estimate of axial and azimuthal induction
        % Calculate thrust coefficient at the streamtube 
        CT = load3Daxial/(0.5*Area*Uinf^2);
        
        % calculate new axial induction, accounting for Glauert's correction
        anew =  ainduction(CT, yaw_angle);
        
        % correct new axial induction with Prandtl's correction
        TSR_input = Omega*Radius/Uinf;
        if UsePrandtl
            [Prandtl, ~, ~] = PrandtlTipRootCorrection(r_R, ...
                rootradius_R, tipradius_R, TSR_input, NBlades, anew);

            if (Prandtl < 0.0001)
                Prandtl = 0.0001; % Avoid divide by zero
            end
        else
            Prandtl = 1;
        end
        % correct estimate of axial induction
        anew = anew/Prandtl; 
        % for improving convergence, weigh current and previous iteration of axial induction
        a = 0.75*a+0.25*anew; 

        % calculate aximuthal induction
        aline = ftan*NBlades/(2*pi*Uinf*(1-a)*Omega*2*(r_R*Radius)^2);
        % correct estimate of azimuthal induction with Prandtl's correction
        aline = aline/Prandtl; 
        
        % Test convergence of solution, by checking convergence of axial induction
        if (abs(a-anew) < Erroriterations) 
            break
        end
    end

    % Calculate Cq
    Cq = 4*aline*(1-a)*Omega*r_R*Radius*r_R;

end