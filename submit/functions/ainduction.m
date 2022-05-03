function a = ainduction(CT, gamma)

    % Calculate the axial induction factor given the thrust coefficient CT
    % including Glauert's correction, and using an approximation for the
    % yawed case

    % Determine an angle after which using the yawed Glauert correction
    % formula
    min_yaw = 5;        % Minimum yaw angle [deg]

    % Preallocate the axial induction factor array
    a = zeros(1,size(CT,2));

    CT1 = 1.816;
    CT2 = 2*sqrt(CT1) - CT1;

    a(CT >= CT2) = 1 + (CT(CT >= CT2) - CT1)/(4*(sqrt(CT1) - 1));
    a(CT < CT2) = 0.5 - 0.5*sqrt(1-CT(CT < CT2));

    % Next we apply an equation given in the KTH-DTU MSc Thesis, it does
    % not make much difference relative to the order of magnitude of the
    % error in the calculation. It can be omitted or checked if something
    % else works better
    if gamma > min_yaw % Go for the yawed case
        ac = 0.2; % From Spera, 1994
        a(a>0.3) = (CT(a>0.3)/4 - ac^2)/(1-2*ac);
    end
end