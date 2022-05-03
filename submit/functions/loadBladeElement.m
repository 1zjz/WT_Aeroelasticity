function [fnorm,ftan,gamma, alpha] = loadBladeElement(vnorm,vtan,chord,...
                              twist, polar_alpha, polar_cl, polar_cd, pitch)                          

%   Calculates the loads in the blade element
%   Inputs:
%   vnorm - normal velocity (Axial velocity) 
%   vtan - tangential velocity
%   chord - chord of the blade element
%   twist - twist angle of the blade element
%   polat_alpha, cl and cd - airfoil coefficients and aoa

    % Magnitude of the square velocity
    vmag2 = vnorm^2+vtan^2;
    % Calculate inflow angle (in rad)
    inflowangle = atand(vnorm/vtan);
    % Add the inflow angle to the twist to get the angle of attack (in deg)
    alpha = - twist + inflowangle + pitch;
    % Use the effective angle of attack to get the lift and drag
    % coefficients cl and cd
    if alpha > polar_alpha(end) 
        % The angle of attack of the blade element is higher than the
        % maximum, then assign the last value of cl, cd
        cl = polar_cl(end); cd = polar_cd(end);
    elseif alpha < polar_alpha(1)
        % The angle of attack of the blade element is smaller than the
        % minimum, then assign the first value of cl, cd
        cl = polar_cl(1); cd = polar_cd(1);
    else
        cl = interp1(polar_alpha,polar_cl,alpha);
        cd = interp1(polar_alpha,polar_cd,alpha);
    end
    % Calculate the lift and drag forces from the cl and cd
    lift = 0.5*vmag2*cl*chord;
    drag = 0.5*vmag2*cd*chord;
    % Decompose lift and drag into normal and tangential components
    fnorm = lift*cosd(inflowangle)+drag*sind(inflowangle);
    ftan =  lift*sind(inflowangle)-drag*cosd(inflowangle);
    % 
    gamma = 0.5*sqrt(vmag2)*cl*chord;

end