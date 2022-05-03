function [CT] = getCT(a, glc, gamma)

    % Calculate the thrust coefficient CT for a given axial induction factor
    % and yaw angle (gamma). The flag 'glc' is used for determining if using
    % the Glauert correction
    
    % The input axial induction factor can be scalar or array
    
    % Hardcoded decision for using yaw case
    min_yaw = 5;        % Minimum yaw angle [deg]

    % Calculate the thrust coefficient (valid for both yaw and non-yawed
    % cases)
    CT = 4*a.*(cosd(gamma) - a);

    if glc == 1 % Use Glauert correction
        if abs(gamma) < min_yaw % Go with the non-yawed case
            CT1 = 1.816;
            a1  = 1 - sqrt(CT1)/2;
            CT(a>a1) = CT1 - 4*(sqrt(CT1) - 1)*(1-a(a>a1));
        else % Go with the yawed case
            % Modify all the values of the thrust coefficient
            CT = 4*a.*sqrt(1 - a.*(2*cosd(gamma) - a));
        end
    end
end