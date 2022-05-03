function [a, a_line, fnorm, ftan, gamma, infangle, aoa_be,clArray, cdArray] =          ...
    solveAzimuthalBEM(N, Nazim, r_blade, chord_blade, twist_blade, ...
    psi_annuli, Uinf, r_R, Root_R, Tip_R ,omega, R, B,             ...
    aoa, cl, cd, yaw_angle, theta, UsePrandtl)

% This functions solves the BEM Method over an azimuthally discretised
% streamtube.
% -------------------------------------------------------------------------
% Inputs
% N           --- Number of blade elements                          [#]
% Nazim       --- Number of azimuthal sectors                       [#]
% r_blade     --- Array containing the blade element discretisation [m]
% chord_blade --- Array containing the chord of each blade element  [m]
% twist_blade --- Array containing the twist of each blade element  [deg]
% psi_annuli  --- Array containing the azimuthal discretisation     [deg]
% Uinf
% r_R
% Root_R
% Tip_R
% omega
% R
% B
% aoa
% cl
% cd
% yaw_angle
% theta
% 
% -------------------------------------------------------------------------
% Intermediate variables
% ii, jj --- Loop variables
% chord --- Chord value of the blade element ii
% twist --- Twist value of the blade element ii
% psi   --- Value of the azimuthal position element jj
% -------------------------------------------------------------------------
% Outputs
% a
% a_line
% fnorm
% ftan
% gamma
% infangle

% *** % *** % *** % *** % *** % *** % *** % *** % *** % *** % *** % *** % *

% Preallocate the variables to store the values in the loops
% Output preallocation
a        = zeros(N-1, 1);         % Mean axial induction factor
a_line   = zeros(N-1, 1);         % Mean tangential induction factor
fnorm    = zeros(N-1, 1);         % Mean normal force per unit length
ftan     = zeros(N-1, 1);         % Mean tangential force per unit length
gamma    = zeros(N-1, 1);         % Mean circulation component
infangle = zeros(N-1, 1);         % Inflow angle
aoa_be   = zeros(N-1, 1);         % Angle of attack of the blade element

% Intermediate values preallocation
a_az        = zeros(N-1, Nazim);     % Axial induction factor for each blade 
                                     % element and azimuthal position
aline_az    = zeros(N-1, Nazim);     % Tangential induction factor
fnormaz     = zeros(N-1, Nazim);     % Normal force component
ftanaz      = zeros(N-1, Nazim);     % Tangential force component
gammaaz     = zeros(N-1, Nazim);     % Circulation component
infangleaz  = zeros(N-1, Nazim);     % Inflow angle
aoa_beaz    = zeros(N-1, Nazim);     % Angle of attack

% The inner loops operate the Blade Element Momentum (BEM) Method 
for ii = 1:N-1            % BLADE ELEMENT LOOP

    % Calculate by interpolation the chord and twist of the blade element
    chord = interp1(r_blade, chord_blade, (r_blade(ii)+r_blade(ii+1))/2);
    twist = interp1(r_blade, twist_blade, (r_blade(ii)+r_blade(ii+1))/2);

    for jj = 1:Nazim    % AZIMUTHAL POSITION LOOP
        % Choose the corresponding azimuthal position
        psi = (psi_annuli(jj+1)+psi_annuli(jj))/2;
        % Call the solve stream tube function
        [a_az(ii,jj), aline_az(ii,jj), fnormaz(ii,jj), ftanaz(ii,jj),  ...
            gammaaz(ii,jj), infangleaz(ii,jj), aoa_beaz(ii,jj), cl_array(ii,jj), cd_array(ii,jj)] =      ...
            solveStreamtube2(Uinf, r_R(ii), r_R(ii+1), Root_R, Tip_R,   ...
            omega, R, B, chord, twist, aoa, cl, cd, yaw_angle, psi,    ...
            theta, UsePrandtl);

    end                 % END AZIMUTHAL POSITION LOOP

    % Calculate the mean values of induction factors and laods
    a(ii)        = mean(a_az(ii,:));
    a_line(ii)   = mean(aline_az(ii,:));
    fnorm(ii)    = mean(fnormaz(ii,:));
    ftan(ii)     = mean(ftanaz(ii,:));
    gamma(ii)    = mean(gammaaz(ii,:));
    infangle(ii) = mean(infangleaz(ii,:));
    aoa_be(ii)   = mean(aoa_beaz(ii,:));
    clArray(ii)  = mean(cl_array(ii,:));
    cdArray(ii)  = mean(cd_array(ii,:));

end                     % END BLADE ELEMENT LOOP

end