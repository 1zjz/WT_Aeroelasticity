%% New BEM for like the millionth and 1 time
clear; close all; clc; 

% Add important folders to the path to call functions and open inputs
functions_folder = ".\functions";
input_folder     = ".\inputs";
addpath(functions_folder);
addpath(input_folder);

% Blade Geometry
B = 3;                % Number of blades                       [#]
R = 50;               % Bladius                                [m]
RootLocation_R = 0.2; % Location of the root over the radius   [-]
TipLocation_R  = 1;   % Location of tip over the radius        [-]
r0 = 0.2*R;           % Radius of first element                [m]
theta = 2;            % Pitch angle                            [deg]
% Initial streamtube discretization
N = 81;               % Number of blade elements               [#]
Nazim = 36;           % Number of azimuthal positions          [#]
% Freestream conditions
rho = 1.225;          % Density of air                         [kg/m3]
Uinf = 10;            % Unperturbed wind speed                 [m/s]
p_inf = 1;            % Freestream pressure                    [Pa]
% Corrections used
glc = 1;              % We want the glauert correction         [-]
UsePrandtl = 1;       % Flag for using Prandtl correction      [-]

% Calculate radial distribution (linearly spaced)
r_blade = linspace(r0, R, N);
% Calculate twist distribution
twist_blade = 14*(1-r_blade/R);
% Calculate chord distribution 
chord_blade = 3*(1-r_blade/R) + 1;



%% Tasks
% Parameters task 1 -------------------------------------------------------
U0 = 10;                    % Downstream wind velocity [m/s]
TSR = 10;                   % range of TSR
file_name = 'DU95W180.csv'; % cl,cd,cm,alpha data for airfoil
omega = U0*TSR/R;           % Rotor speeds for TSR

% Parameters task 2 -------------------------------------------------------
yaw_angle = 0;                % Yaw angle [deg] 
% tsr_2 = 8;                  % Tip-speed ratio for task 2 [-]
% 
% % Parameters task 3 -------------------------------------------------------
% ct_3 = 0.75;                % Limiting Thrust Coefficient 
% tsr_3 = 8;                  % Design-tip-speed ratio
% area = 2*pi*R;
% 
% % Plot flag
% plotChecks = 0;

%% Load Polar Curves of the aerofoil and plot them
polar_data = table2array(readtable(file_name)); 
aoa = polar_data(:,1); %[deg]
cl = polar_data(:,2);
cd = polar_data(:,3);
cm = polar_data(:,4);

% Calculate maximum lift-to-drag ratio
Cl_Cd = cl./cd; [ClCdmax, iiCldesign] = max(Cl_Cd);
Cldesign = cl(iiCldesign); aoa_design = aoa(iiCldesign);


%% Plot the CT curve for non-yawed and yawed case
% Using the axial induction factor and the Glauert correction
% Series of induction factor
a_az = linspace(-0.5, 1, 100);

for ii = 1:length(yaw_angle)
    CTmom     = getCT(a_az, 0, yaw_angle(ii)); % CT without correction
    CTglauert = getCT(a_az, 1, yaw_angle(ii)); % CT with Glauert's correction

    % Calculate the induction factor based on the thrust coefficient
    a_back = ainduction(CTglauert, yaw_angle(ii));
    % Check the difference with the actual axial induction factor
    error_a = a_az - a_back;

end

%% Prandtl Tip Correction Function Test
% Axial induction factor to test the Prandtl's tip correction
a_test = 0.5;
% Define an array of radial positions, from root to tip
r_R = linspace(RootLocation_R, TipLocation_R, N);
% Define the centroides of the radial elements
r_Rc = 0.5*(r_R(2:end) + r_R(1:end-1));
% Call Prandtl Tip Root Correction function
[Prandtl, Prandtltip, Prandtlroot] = ...
    PrandtlTipRootCorrection(r_R, RootLocation_R, TipLocation_R, TSR(1),...
    B, a_test);

%% Call BEM Algorithm
% Solve BEM model (Iteratively go through the tip-speed ratios, blade
% elements and azimuthal positions)

% Define the azimuthal position array 
psi_annuli = linspace(0, 360, Nazim+1); % Azimuthal position [deg]
% Define the non-dimensional array position
clear r_R
r_R = r_blade / R;
% Discretization step in radial direction
dr_R = r_R(2:end) - r_R(1:end-1);
dr   = dr_R * R;
% Area annuli elements 
% dA   = (r_R(2:end).^2 - r_R(1:end-1).^2)*pi*R^2; % ----> NOT USED

%% Task 1 · BEM for different TSR of a Wind Turbine

% Preallocate the variables to store them in the loop
CT       = zeros(length(TSR), 1);   % Thrust coefficient for each TSR
CT_alt   = zeros(length(TSR), 1);   % Thrust coefficient for each TSR
CP       = zeros(length(TSR), 1);   % Power coefficient for each TSR
CQ       = zeros(length(TSR), 1);   % Torque coefficient for each TSR
CQ_alt   = zeros(length(TSR), 1);   % Torque coefficient for each TSR
Thrust   = zeros(N-1, length(TSR)); % Thrust for each TSR
Torque   = zeros(N-1, length(TSR)); % Torque for each TSR
Power    = zeros(N-1, length(TSR)); % Power for each TSR
Q        = zeros(length(TSR), 1);   % Torque for each TSR
P        = zeros(length(TSR), 1);   % Power for each TSR
T        = zeros(length(TSR), 1);   % Thrust for each TSR
a        = zeros(N-1, length(TSR)); % Axial induction factor
aline    = zeros(N-1, length(TSR)); % Tangential induction factor
fnorm    = zeros(N-1, length(TSR)); % Normal force component
ftan     = zeros(N-1, length(TSR)); % Tangential force component
gamma    = zeros(N-1, length(TSR)); % Circulation 
infang   = zeros(N-1, length(TSR)); % Inflow angle
aoa_be   = zeros(N-1, length(TSR)); % Angle of attack of the blade element
Cq       = zeros(N-1, length(TSR)); % Torque coeff of the blade element

for ii = 1:length(TSR)      % TIP-SPEED RATIO LOOP

    % Call the BEM code for azimuthal discretisation solution
    [a(:,ii), aline(:,ii), fnorm(:,ii), ftan(:,ii), ...
        gamma(:,ii), infang(:,ii), aoa_be(:,ii), Cq(:,ii)] = ...
    solveAzimuthalBEM(N, Nazim, r_blade, chord_blade, twist_blade, ...
    psi_annuli, Uinf, r_R, RootLocation_R, TipLocation_R ,omega(ii), R, B, ...
    aoa, cl, cd, yaw_angle(1), theta, UsePrandtl);


    % Calculate the overall power and thrust coefficients
    CT(ii) = sum(dr.*fnorm(:,ii)'*B/(0.5*Uinf^2*pi*R^2));
    CQ(ii) = sum(dr.*ftan(:,ii)'.*r_Rc*R*B/(0.5*Uinf^2*pi*R^3));
    CP(ii) = sum(dr.*ftan(:,ii)'.*r_Rc*R*B*omega(ii)/(0.5*Uinf^3*pi*R^2));
    Cq(:,ii) = dr.*ftan(:,ii)'.*r_Rc*R*B/(0.5*Uinf^2*pi*R^3);
    
    % Thrust, torque and power distribution
    Thrust(:,ii) = fnorm(:,ii)'*B;
    Torque(:,ii) = ftan(:,ii)'.*r_Rc*R*B;
    Power(:,ii) = ftan(:,ii)'.*r_Rc*R*B*omega(ii);
    % Total power, thrust and torque
    P(ii)  = CP(ii) * (0.5*Uinf^3*pi*R^2);
    Q(ii)  = P(ii) / omega(ii);
    T(ii)  = CT(ii) * 0.5*Uinf^2*pi*R^2;

end                         % END TIP-SPEED RATIO LOOP


% %% Stagnation Enthalpy
% 
% A = zeros((N),length(TSR));
% A1 = zeros((N),length(TSR));
% A4 = zeros((N),length(TSR));
% CT_r = zeros((N),length(TSR));
% h1 = zeros((N-1),length(TSR));
% h2 = zeros((N-1),length(TSR));
% h3 = zeros((N-1),length(TSR));
% h4 = zeros((N-1),length(TSR));
% r4 = zeros(N,length(TSR));
% r1 = zeros(N,length(TSR));
% r4_c = zeros(N-1,length(TSR));
% r1_c = zeros(N-1,length(TSR));
% r_c = zeros(N-1,1);
% 
% count = zeros(1, length(TSR));
% 
% for jj = 1:length(TSR)
%     count(jj) = 1;
%    
%     %Velcity at the root at actuator disk
%     U_root_cyl = Uinf*(1-a(2,jj));
%     %Velocity at the root of the wake
%     %assuming a at root is same as a at first radial position on blade
%     U_root_wake = Uinf*(1-2*a(2,jj));
%     %Velocity at root of downstream
%     U_root_down = Uinf;
% 
%     %Area of the root at actuator disk
%     A_root_cyl = pi*r_blade(3)^2; 
%     %Area of the root of the wake
%     A_root_wake = U_root_cyl/U_root_wake*A_root_cyl;
%     %Area of root upstrea
%     A_root_up = U_root_cyl/U_root_down*A_root_cyl;
%     
%     %First radial position at wake
%     r4(1,jj) = sqrt(A_root_wake/pi);
%     %First radial position upstream
%     r1(1,jj) = sqrt(A_root_up/pi);
% 
% 
% 
%     for ii = 3:N-1
% 
%         if a(ii,jj) < 0.5
%             count(jj) = count(jj) +1;
%         %Calculation of Stagnation enthalpy
%         % Area of the streamtube element at the actuator disk
%         A(ii,jj) = pi*((r_blade(ii+1))^2 - (r_blade(ii))^2); 
%         % CT for at each radial position on the actuator disk
%         CT_r(ii,jj) = dr(ii).*fnorm(ii,jj)'*B/(0.5*Uinf^2*pi*R^2);
%         % Calculate stagnation enthalpy
%         [h1(ii-2,jj), h2(ii-2,jj), h3(ii-2,jj), h4(ii-2,jj)]...
%             = stagnationEnthalpy(Uinf, rho, p_inf, a(ii,jj), CT_r(ii,jj), A(ii,jj));
%         
%         
%         % Calculation of radial position and streamtube expansion
%         % Mean radial position
%         r_c(ii) = (r_blade(ii) + r_blade(ii+1))/2;
% 
%         % Area at streamtube inlet 
%         A1(ii-2,jj) = (1-a(ii,jj))*A(ii,jj); 
%         r1(ii-1,jj) = sqrt((A1(ii-2,jj)/pi) + r1(ii-2,jj)^2); % + r1(ii-2,jj);
%         r1_c(ii-2,jj) = (r1(ii-1,jj) + r1(ii-2,jj))/2;
%         % Area at streamtube outlet 
%         A4(ii-2,jj) = ((1-a(ii,jj))/(1-(2*a(ii,jj))))*A(ii,jj);
%         r4(ii-1,jj) = sqrt((A4(ii-2,jj)./pi) + r4(ii-2,jj)^2); % + r4(ii-2,jj);
%         r4_c(ii-2,jj) = (r4(ii-1,jj) + r4(ii-2,jj))/2;
% 
%         else 
%             break
%         end
%     end
%     
% end

% %% Task 2 · BEM for different yaw angles of a Wind Turbine
% 
% % Number of yaw angles to analyse
% Nyaw = length(yaw_angle);
% 
% % Preallocate the variables to store them in the loop
% CT_yaw    = zeros(Nyaw, 1); % Thrust coefficient for each TSR
% CP_yaw    = zeros(Nyaw, 1); % Power coefficient for each TSR
% 
% a        = zeros(N-1, Nyaw); % Axial induction factor
% aline    = zeros(N-1, Nyaw); % Tangential induction factor
% fnorm    = zeros(N-1, Nyaw); % Normal force component
% ftan     = zeros(N-1, Nyaw); % Tangential force component
% gamma    = zeros(N-1, Nyaw); % Circulation 
% infang   = zeros(N-1, Nyaw); % Inflow angle
% aoa_be   = zeros(N-1, Nyaw); % Angle of attack
% 
% % Fix the tip speed ratio and rotational speed
% clear TSR; TSR = 8; 
% clear omega; omega = TSR * Uinf / R;
% 
% for ii = 1:length(yaw_angle) % YAW ANGLE LOOP
% 
% % Call the BEM code for azimuthal discretisation solution
%     [a(:,ii), aline(:,ii), fnorm(:,ii), ftan(:,ii), ...
%         gamma(:,ii), infang(:,ii), aoa_be(:,ii)] = ...
%     solveAzimuthalBEM(N, Nazim, r_blade, chord_blade, twist_blade, ...
%     psi_annuli, Uinf, r_R, RootLocation_R, TipLocation_R ,omega, R, B, ...
%     aoa, cl, cd, yaw_angle(ii), theta, UsePrandtl);
% 
%     % Plot the mean axial and tangential induction factors along the blade
%     figure("Name",strcat('Induction factors for yaw ', ...
%             num2str(yaw_angle(ii)), ' and TSR ', num2str(TSR)));
%     hold on
%     plot(r_blade(1:end-1), a, 'DisplayName','$a$');
%     plot(r_blade(1:end-1), aline, 'DisplayName','$a^{\prime}$');
%     grid on
%     xlabel('Radial $r$ [m]','Interpreter','latex')
%     ylabel('Induction factor [-]','Interpreter','latex')
%     legend('Location','best','Interpreter','latex')
% 
%     % Plot the mean normal and tangential forces along the blade
%     figure("Name",strcat('Mean loads for yaw ', ...
%             num2str(yaw_angle(ii)), ' and TSR ', num2str(TSR)))
%     hold on
%     plot(r_blade(1:end-1), fnorm, 'DisplayName','$f_{norm}$');
%     plot(r_blade(1:end-1), ftan, 'DisplayName','$f_{tan}$');
%     grid on
%     xlabel('Radial position $r$ [m]','Interpreter','latex')
%     ylabel('Load [N/m]','Interpreter','latex')
%     legend('Location','best','Interpreter','latex')
%     
%     % ##############################################################################################
%     % TEST FNORM AND FTAN AS SUBPLOTS
%     figure(111)
%     title(strcat('TEST: Mean loads for yaw ', ...
%            num2str(yaw_angle(ii)), ' and TSR ', num2str(TSR)))
%     hold on
%     subplot(2,1,1)
%     plot(r_blade(1:end-1), fnorm(:,ii), 'DisplayName',['$f_{norm}$ at yaw ', ...
%            num2str(yaw_angle(ii)), ' and TSR ', num2str(TSR)']);
%     hold on
%     grid on
%     xlabel('Radial position $r$ [m]','Interpreter','latex')
%     ylabel('Load [N/m]','Interpreter','latex')
%     legend('Location','best','Interpreter','latex')
%     subplot(2,1,2)
%     plot(r_blade(1:end-1), ftan(:,ii), 'DisplayName',['$f_{tan}$ at yaw ', ...
%            num2str(yaw_angle(ii)), ' and TSR ', num2str(TSR)']);
%     hold on
%     grid on
%     xlabel('Radial position $r$ [m]','Interpreter','latex')
%     ylabel('Load [N/m]','Interpreter','latex')
%     legend('Location','best','Interpreter','latex')
%     
%     
%     % #############################################################################################
% 
%     % Calculate the overall power and thrust coefficients
%     CT_yaw(ii) = sum(dr.*fnorm(:,ii)'*B/(0.5*Uinf^2*pi*R^2));
%     CP_yaw(ii) = sum(dr.*ftan(:,ii)'.*r_Rc*B*R*omega/(0.5*Uinf^3*pi*R^2));
% 
% end % END YAW ANGLE LOOP
% 
% % Centered azimuthal positions
% psi_c = (psi_annuli(1:end-1) + psi_annuli(2:end))*0.5; 

% %% Analysis of the use of tip and root correction
% % Prandtl flag
% UsePrandtl = [1 0];
% 
% % Preallocate variables
% CP_PrandtlComp = zeros(2, 1); % Power coefficient
% CT_PrandtlComp = zeros(2, 1); % Thrust coefficient
% a              = zeros(N-1, 2); % Axial induction factor
% aline          = zeros(N-1, 2); % Tangential induction factor
% fnorm          = zeros(N-1, 2); % Normal force component
% ftan           = zeros(N-1, 2); % Tangential force component
% gamma          = zeros(N-1, 2); % Circulation 
% infang         = zeros(N-1, 2); % Inflow angle
% aoa_be         = zeros(N-1, 2); % Angle of attack
% 
% for ii = 1:2
% 
%     % Call the BEM code for azimuthal discretisation solution
%     [a(:,ii), aline(:,ii), fnorm(:,ii), ftan(:,ii), ...
%         gamma(:,ii), infang(:,ii), aoa_be(:,ii)] = ...
%     solveAzimuthalBEM(N, Nazim, r_blade, chord_blade, twist_blade, ...
%     psi_annuli, Uinf, r_R, RootLocation_R, TipLocation_R ,omega, R, B, ...
%     aoa, cl, cd, yaw_angle(1), theta, UsePrandtl(ii));
% 
% end