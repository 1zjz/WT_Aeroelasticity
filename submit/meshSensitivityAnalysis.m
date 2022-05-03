%% Sensitivity analysis of the BEM Code

% This script deals with the discretisation of the blade. It compares the
% use of more finely discretised stream tubes in the radial and azimuthal
% directions. Also it compares the use of differently spaced nodes to deal
% with the high gradients at the root and tip of the blade

clear; close all; clc; 

%% Parameters definition
% Add important folders to the path to call functions and open inputs
functions_folder = ".\functions";
input_folder     = ".\inputs";
addpath(functions_folder);
addpath(input_folder);

B = 3;                % Number of blades                       [#]
R = 50;               % Bladius                                [m]
Nazim = 36;           % Number of azimuthal positions          [#]
rho = 1.225;          % Density of air                         [kg/m3]
RootLocation_R = 0.2; % Location of the root over the radius   [-]
TipLocation_R  = 1;   % Location of tip over the radius        [-]
r0 = 0.2*R;           % Radius of first element                [m]
theta = 2;            % Pitch angle                            [deg]
glc = 1;              % We want the glauert correction         [-]
Uinf = 1;             % Unperturbed wind speed                 [m/s]
TSR  = 8;             % Tip-Speed Ratio                        [-]
yaw_angle = 0;        % Yaw angle                              [deg]
UsePrandtl = 1;       % Flag to use Prandtl's tip and root correction

% Calculate rotor speed for the current case
omega = TSR * Uinf / R;
% Load polar data
file_name = 'DU95W180.csv'; % cl,cd,cm,alpha data for airfoil
polar_data = table2array(readtable(file_name)); 
aoa = polar_data(:,1); %[deg]
cl = polar_data(:,2);
cd = polar_data(:,3);
cm = polar_data(:,4);

% Define the azimuthal position array 
psi_annuli = linspace(0, 360, Nazim+1); % Azimuthal position [deg]

%% Influence of number of annuli (for equispaced case)

StepNodes = 5; MinNodes = 20; MaxNodes = 120;
Nradii = MinNodes:StepNodes:MaxNodes;
N_Nnodes = length(Nradii);

% Preallocate the variables to store them in the loop
CT_linR  = zeros(1, N_Nnodes);   % Thrust coefficient for each TSR
CP_linR  = zeros(1, N_Nnodes);   % Power coefficient for each TSR
% Preallocate intermediate values
a        = zeros(MaxNodes, 1); % Axial induction factor
aline    = zeros(MaxNodes, 1); % Tangential induction factor
fnorm    = zeros(1, MaxNodes); % Normal force component
ftan     = zeros(1, MaxNodes); % Tangential force component
gamma    = zeros(MaxNodes, 1); % Circulation 
infang   = zeros(MaxNodes, 1); % Inflow angle
aoa_be   = zeros(MaxNodes, 1); % Angle of attack of the blade element
time_ls  = zeros(N_Nnodes, 1);
Ntimer   = 5;

for ii = 1:N_Nnodes

    % Number of blade elements for the current case
    Nbe = Nradii(ii);
    
    tic
    for jj = 1:Ntimer % Loop for timing
    % Calculate radial distribution (linearly spaced)
    r_blade = linspace(r0, R, Nbe);
    % Calculate twist distribution
    twist_blade = 14*(1-r_blade/R);
    % Calculate chord distribution
    chord_blade = 3*(1-r_blade/R) + 1;
    % Define the non-dimensional array position
    r_R = r_blade / R;
    % Define the centroides of the radial elements
    r_Rc = 0.5*(r_R(2:end) + r_R(1:end-1));
    % Discretization step in radial direction
    dr_R = r_R(2:end) - r_R(1:end-1);
    dr   = dr_R * R;

    % Call the BEM code for azimuthal discretisation solution
    [~, ~, fnorm(1:Nbe-1), ftan(1:Nbe-1), ~, ~, ~] = ...
    solveAzimuthalBEM(Nbe, Nazim, r_blade, chord_blade, twist_blade, ...
    psi_annuli, Uinf, r_R, RootLocation_R, TipLocation_R ,omega, R, B, ...
    aoa, cl, cd, yaw_angle, theta, UsePrandtl);

    % Calculate the overall power and thrust coefficients
    CT_linR(ii) = sum(dr.*fnorm(1:Nbe-1)*B/(0.5*Uinf^2*pi*R^2));
    CP_linR(ii) = sum(dr.*ftan(1:Nbe-1).*r_Rc*B*R*omega/(0.5*Uinf^3*pi*R^2));
    
    end 
    time_ls(ii) = toc/Ntimer;
    % Clear variables to not superimpose arrays of different size
    clear a aline fnorm ftan gamma infang aoa_be r_blade twist_blade ...
        chord_blade r_R dr r_Rc dr_R
end

errCP_linR = (CP_linR - CP_linR(end))./CP_linR(end);
errCT_linR = (CT_linR - CT_linR(end))./CT_linR(end);

%% Plot the convergence of CP and CT
figure("Name","Convergence of CP")
plot(Nradii, CP_linR, 'k-o');
xlabel('Number of blade elements, $N$ [-]', 'Interpreter','latex')
ylabel('Power Coefficient, $C_{P}$ [-]', 'Interpreter','latex')
grid on

figure("Name","Error of CP")
yyaxis left
    plot(Nradii, errCP_linR*100, '-o');
    ylabel('Power Coefficient Error $\varepsilon_{C_{P}}$ [\%]', 'Interpreter','latex')
yyaxis right
    plot(Nradii, time_ls, '-o');
    ylabel('Computation time $t$ [s]', 'Interpreter','latex')
xlabel('Number of blade elements, $N$ [-]', 'Interpreter','latex')

grid on

figure("Name","Convergence of CT")
plot(Nradii, CT_linR, 'k-o');
xlabel('Number of blade elements, $N$ [-]', 'Interpreter','latex')
ylabel('Thrust Coefficient, $C_{T}$ [-]', 'Interpreter','latex')
grid on

figure("Name","Error of CT")
yyaxis left
    plot(Nradii, errCP_linR*100, '-o');
    ylabel('Thrust Coefficient Error $\varepsilon_{C_{T}}$ [\%]', ...
        'Interpreter','latex')
yyaxis right
    plot(Nradii, time_ls, '-o');
    ylabel('Computation time $t$ [s]', 'Interpreter','latex')
xlabel('Number of blade elements, $N$ [-]', 'Interpreter','latex')
grid on

%% Influence of number of azimuthal positions for a yaw case

% Set a different yaw angle
yaw_angle = 15;
% Set minimum and maximum values of azimuthal discretisation
MinAzim = 20; MaxAzim = 150; StepAzim = 20;
% Obtain a linear distribution of discretization and the number of meshes
AzValues = MinAzim:StepAzim:MaxAzim;
NAz = length(AzValues);

% Set the number of blade elements
Nbe = 100;
% Calculate radial distribution (linearly spaced)
r_blade = linspace(r0, R, Nbe);
% Calculate twist distribution
twist_blade = 14*(1-r_blade/R);
% Calculate chord distribution
chord_blade = 3*(1-r_blade/R) + 1;
% Define the non-dimensional array position
r_R = r_blade / R;
% Define the centroides of the radial elements
r_Rc = 0.5*(r_R(2:end) + r_R(1:end-1));
% Discretization step in radial direction
dr_R = r_R(2:end) - r_R(1:end-1);
dr   = dr_R * R;

% Preallocate the variables to store them in the loop
CT       = zeros(1, NAz);   % Thrust coefficient for each TSR
CP       = zeros(1, NAz);   % Power coefficient for each TSR
% Preallocate intermediate values
fnorm    = zeros(1, Nbe-1); % Normal force component
ftan     = zeros(1, Nbe-1); % Tangential force component
time_az  = zeros(1, NAz);
Ntimer   = 5;

for ii = 1:NAz

    tic
    for jj = 1:Ntimer
    % Define the azimuthal position array 
    psi_annuli = linspace(0, 360, AzValues(ii)+1); % Azimuthal position [deg]

    % Call the BEM code for azimuthal discretisation solution
    [~, ~, fnorm, ftan, ~, ~, ~] = ...
    solveAzimuthalBEM(Nbe, AzValues(ii), r_blade, chord_blade, twist_blade, ...
    psi_annuli, Uinf, r_R, RootLocation_R, TipLocation_R ,omega, R, B, ...
    aoa, cl, cd, yaw_angle, theta, UsePrandtl);

    % Calculate the overall power and thrust coefficients
    CT(ii) = sum(dr.*fnorm'*B/(0.5*Uinf^2*pi*R^2));
    CP(ii) = sum(dr.*ftan'.*r_Rc*B*R*omega/(0.5*Uinf^3*pi*R^2));

    end
    time_az(ii) = toc / Ntimer;

end

%% Plots for azimuthal discretisation

errCP_azim = (CP - CP(end))./CP(end);
errCT_azim = (CT - CT(end))./CT(end);

% Plot the convergence of CP and CT
figure("Name","Convergence of CP")
plot(AzValues, CP, 'k-o');
xlabel('Number of azimuthal elements, $N_{az}$ [-]', 'Interpreter','latex')
ylabel('Power Coefficient, $C_{P}$ [-]', 'Interpreter','latex')
grid on

figure("Name","Error of CP")
yyaxis left
    plot(AzValues, errCP_azim*100, '-o');
    ylabel('Power Coefficient Error $\varepsilon_{C_{P}}$ [\%]', 'Interpreter','latex')
yyaxis right
    plot(AzValues, time_az, '-o');
    ylabel('Computation time $t$ [s]', 'Interpreter','latex')
xlabel('Number of azimuthal elements, $N_{az}$ [-]', 'Interpreter','latex')

grid on

figure("Name","Convergence of CT")
plot(AzValues, CT, 'k-o');
xlabel('Number of azimuthal elements, $N_{az}$ [-]', 'Interpreter','latex')
ylabel('Thrust Coefficient, $C_{T}$ [-]', 'Interpreter','latex')
grid on

figure("Name","Error of CT")
yyaxis left
    plot(AzValues, errCP_azim*100, '-o');
    ylabel('Thrust Coefficient Error $\varepsilon_{C_{T}}$ [\%]', ...
        'Interpreter','latex')
yyaxis right
    plot(AzValues, time_az, '-o');
    ylabel('Computation time $t$ [s]', 'Interpreter','latex')
xlabel('Number of azimuthal elements, $N_{az}$ [-]', 'Interpreter','latex')
grid on

%% Influence of discretisation kind (linearly spaced of Chebyshev) 

% Set the yaw angle and the azimuthal discretization
Nazim = 150; yaw_angle = 0;
psi_annuli = linspace(0, 360, Nazim+1); % Azimuthal position [deg]

% Preallocate CP and CT arrays
CP_cosine = zeros(N_Nnodes,1);
CT_cosine = zeros(N_Nnodes,1);
time_cos  = zeros(N_Nnodes,1);


for ii = 1:N_Nnodes
    
    tic
    for jj = 1:Ntimer
    % Using a cosine method for non-linear discretization of the radial
    % direction
    % Set the number of blade elements
    Nbe = Nradii(ii);
    % Calculate radial distribution (linearly spaced)
    points = linspace(-pi/2, 0, Nbe);
    r_blade = RootLocation_R*R + (R - RootLocation_R*R) * cos(points);
    % Calculate twist distribution
    twist_blade = 14*(1-r_blade/R);
    % Calculate chord distribution
    chord_blade = 3*(1-r_blade/R) + 1;
    % Define the non-dimensional array position
    r_R = r_blade / R;
    % Define the centroides of the radial elements
    r_Rc = 0.5*(r_R(2:end) + r_R(1:end-1));
    % Discretization step in radial direction
    dr_R = r_R(2:end) - r_R(1:end-1);
    dr   = dr_R * R;


    % Call the BEM code for azimuthal discretisation solution
    [~, ~, fnorm, ftan, ~, ~, ~] = ...
        solveAzimuthalBEM(Nbe, Nazim, r_blade, chord_blade, twist_blade, ...
        psi_annuli, Uinf, r_R, RootLocation_R, TipLocation_R ,omega, R, B, ...
        aoa, cl, cd, yaw_angle, theta, UsePrandtl);

    % Calculate the overall power and thrust coefficients
    CT_cosine(ii) = sum(dr.*fnorm'*B/(0.5*Uinf^2*pi*R^2));
    CP_cosine(ii) = sum(dr.*ftan'.*r_Rc*B*R*omega/(0.5*Uinf^3*pi*R^2));

    % Clear variables to not superimpose arrays of different size
    clear a aline fnorm ftan gamma infang aoa_be r_blade twist_blade ...
        chord_blade r_R dr r_Rc dr_R points
    end
    time_cos(ii) = toc/Ntimer;
end

%% Plot cosine versus linearly spaced
errCP_cos = (CP_cosine - CP_cosine(end))./CP_cosine(end);
errCT_cos = (CT_cosine - CT_cosine(end))./CT_cosine(end);

% Compare with the linear case
figure("Name","Comparison of discretisation scheme using CP")
hold on
plot(Nradii, CP_cosine, 'k-o', 'DisplayName', 'Cosine')
plot(Nradii, CP_linR, 'k--o', 'DisplayName', 'Linear')
grid on
xlabel('Number of blade elements, $N$ [-]', 'Interpreter','latex')
ylabel('Power Coefficient, $C_{P}$ [-]', 'Interpreter','latex')
legend('Interpreter','latex','Location','best')

figure("Name","Error of CP")
yyaxis left
    plot(Nradii, errCP_cos*100, '-o', 'DisplayName', 'Cosine');
    hold on
    plot(Nradii, errCP_linR*100, '--o', 'DisplayName', 'Linearly-spaced');
    ylabel('Power Coefficient Error $\varepsilon_{C_{P}}$ [\%]', ...
        'Interpreter','latex')
yyaxis right
    plot(Nradii, time_cos, '-o', 'DisplayName', 'Cosine');
    hold on;
    plot(Nradii, time_ls, '--o', 'DisplayName', 'Linearly-spaced');
    ylabel('Computation time $t$ [s]', 'Interpreter','latex')

xlabel('Number of radial elements, $N$ [-]', 'Interpreter','latex')
legend('Interpreter','latex','Location','best')
grid on

figure("Name","Comparison of discretisation scheme using CT")
hold on
plot(Nradii, CT_cosine, 'k-o', 'DisplayName', 'Cosine')
plot(Nradii, CT_linR, 'k--o', 'DisplayName', 'Linear')
grid on
xlabel('Number of blade elements, $N$ [-]', 'Interpreter','latex')
ylabel('Thrust Coefficient, $C_{T}$ [-]', 'Interpreter','latex')
legend('Interpreter','latex','Location','best')

figure("Name","Error of CT")
yyaxis left
    hold on
    plot(Nradii, errCT_cos*100, '-o', 'DisplayName', 'Cosine');
    plot(Nradii, errCT_linR*100, '--o', 'DisplayName', 'Linearly-spaced');
    ylabel('Power Coefficient Error $\varepsilon_{C_{P}}$ [\%]', 'Interpreter','latex')
yyaxis right
    plot(Nradii, time_cos, '-o', 'DisplayName', 'Cosine');
    hold on;
    plot(Nradii, time_ls, '--o', 'DisplayName', 'Linearly-spaced');
    ylabel('Computation time $t$ [s]', 'Interpreter','latex')
xlabel('Number of radial elements, $N$ [-]', 'Interpreter','latex')
legend('Interpreter','latex','Location','best')
grid on
