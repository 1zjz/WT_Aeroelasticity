%% New BEM for like the millionth and 1 time
clear; close all; clc; 

% Add important folders to the path to call functions and open inputs
functions_folder = ".\functions";
input_folder     = ".\inputs";
addpath(functions_folder);
addpath(input_folder);

B = 3;                % Number of blades                       [#]
R = 50;               % Bladius                                [m]
N = 81;               % Number of blade elements               [#]
Nazim = 36;           % Number of azimuthal positions          [#]
rho = 1.225;          % Density of air                         [kg/m3]
RootLocation_R = 0.2; % Location of the root over the radius   [-]
TipLocation_R  = 1;   % Location of tip over the radius        [-]
r0 = 0.2*R;           % Radius of first element                [m]
theta = 2;            % Pitch angle                            [deg]
glc = 1;              % We want the glauert correction         [-]
Uinf = 1;             % Unperturbed wind speed                 [m/s]
p_inf = 1;            % Freestream pressure                    [Pa]
UsePrandtl = 1;       % Flag for using Prandtl correction
% Calculate radial distribution (linearly spaced)
r_blade = linspace(r0, R, N);
% Calculate twist distribution
twist_blade = 14*(1-r_blade/R);
% Calculate chord distribution 
chord_blade = 3*(1-r_blade/R) + 1;

% Parameters task 1
U0 = 1;                     % Downstream wind velocity [m/s]
TSR = 8;             % [6 8 10]; % range of TSR
file_name = 'DU95W180.csv'; % cl,cd,cm,alpha data for airfoil
omega = U0*TSR/R;           % Rotor speeds for TSR

% Parameters task 2 
yaw_angle = 0;      % Yaw angle [deg] 


%% Load Polar Curves of the aerofoil and plot them
polar_data = table2array(readtable(file_name)); 
aoa = polar_data(:,1); %[deg]
cl = polar_data(:,2);
cd = polar_data(:,3);
cm = polar_data(:,4);

% figure('Name', 'Cl curve')
% plot(aoa,cl)
% grid on
% xlabel('Angle of attack $\alpha$ [deg]', 'Interpreter','latex');
% ylabel('Lift Coefficient $C_{l}$ [-]', 'Interpreter','latex');
% 
% figure('Name', 'Cd curve')
% plot(aoa,cd)
% grid on
% xlabel('Angle of attack $\alpha$ [deg]', 'Interpreter','latex');
% ylabel('Drag Coefficient $C_{l}$ [-]', 'Interpreter','latex');
% 
% figure('Name', 'Cm curve')
% plot(aoa,cm)
% grid on
% xlabel('Angle of attack $\alpha$ [deg]', 'Interpreter','latex');
% ylabel('Moment Coefficient $C_{l}$ [-]', 'Interpreter','latex');

% Interesting to add cl vs cd plot


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
CP       = zeros(length(TSR), 1);   % Power coefficient for each TSR
a        = zeros(N-1, length(TSR)); % Axial induction factor
aline    = zeros(N-1, length(TSR)); % Tangential induction factor
fnorm    = zeros(N-1, length(TSR)); % Normal force component
ftan     = zeros(N-1, length(TSR)); % Tangential force component
gamma    = zeros(N-1, length(TSR)); % Circulation 
infang   = zeros(N-1, length(TSR)); % Inflow angle
aoa_be   = zeros(N-1, length(TSR)); % Angle of attack of the blade element

for ii = 1:length(TSR)      % TIP-SPEED RATIO LOOP

    % Call the BEM code for azimuthal discretisation solution
    [a(:,ii), aline(:,ii), fnorm(:,ii), ftan(:,ii), ...
        gamma(:,ii), infang(:,ii), aoa_be(:,ii),clArray(:,ii),cdArray(:,ii)] = ...
    solveAzimuthalBEM2(N, Nazim, r_blade, chord_blade, twist_blade, ...
    psi_annuli, Uinf, r_R, RootLocation_R, TipLocation_R ,omega(ii), R, B, ...
    aoa, cl, cd, yaw_angle(1), theta, UsePrandtl);

    % Calculate the overall power and thrust coefficients
%     CT(ii) = sum(dr.*fnorm(:,ii)'*B/(0.5*Uinf^2*pi*R^2));
%     CP(ii) = sum(dr.*ftan(:,ii)'.*r_Rc*B*R*omega(ii)/(0.5*Uinf^3*pi*R^2));


end                         % END TIP-SPEED RATIO LOOP


%%
cl_cd = clArray ./ cdArray;
OperationalPointArray = [clArray,cdArray];

save('OperationalPoint.mat','OperationalPointArray')


%% Maximized Blade

load('optBlade.mat');


% Calculate chord distribution 
chord_blade = optBlade(1,:);

% Calculate twist distribution
twist_blade = optBlade(2,:);

% Call BEM Algorithm
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

% Task 1 · BEM for different TSR of a Wind Turbine

% Preallocate the variables to store them in the loop
CT       = zeros(length(TSR), 1);   % Thrust coefficient for each TSR
CP       = zeros(length(TSR), 1);   % Power coefficient for each TSR
a        = zeros(N-1, length(TSR)); % Axial induction factor
aline    = zeros(N-1, length(TSR)); % Tangential induction factor
fnorm    = zeros(N-1, length(TSR)); % Normal force component
ftan     = zeros(N-1, length(TSR)); % Tangential force component
gamma    = zeros(N-1, length(TSR)); % Circulation 
infang   = zeros(N-1, length(TSR)); % Inflow angle
aoa_be2   = zeros(N-1, length(TSR)); % Angle of attack of the blade element

for ii = 1:length(TSR)      % TIP-SPEED RATIO LOOP

    % Call the BEM code for azimuthal discretisation solution
    [a(:,ii), aline(:,ii), fnorm(:,ii), ftan(:,ii), ...
        gamma(:,ii), infang(:,ii), aoa_be2(:,ii),clArray(:,ii),cdArray(:,ii)] = ...
    solveAzimuthalBEM2(N, Nazim, r_blade, chord_blade, twist_blade, ...
    psi_annuli, Uinf, r_R, RootLocation_R, TipLocation_R ,omega(ii), R, B, ...
    aoa, cl, cd, yaw_angle(1), theta, UsePrandtl);

    % Calculate the overall power and thrust coefficients
%     CT(ii) = sum(dr.*fnorm(:,ii)'*B/(0.5*Uinf^2*pi*R^2));
%     CP(ii) = sum(dr.*ftan(:,ii)'.*r_Rc*B*R*omega(ii)/(0.5*Uinf^3*pi*R^2));


end                         % END TIP-SPEED RATIO LOOP


OperationalPointArray2 = [clArray,cdArray];

save('OperationalPoint2.mat','OperationalPointArray2')


%% Load Files for Blades

load("OperationalPoint.mat");
load("OperationalPoint2.mat");
load("optBlade.mat");
load("optAngles.mat");



orig_cl     = OperationalPointArray(:,1);
orig_cd     = OperationalPointArray(:,2);
orig_clcd   = orig_cl ./ orig_cd;

opt_cl      = OperationalPointArray2(:,1);
opt_cd      = OperationalPointArray2(:,2);
opt_clcd    = opt_cl ./ opt_cd;


r_blade = linspace(r0, R, N);

orig_twist  = 14*(1-r_blade/R);
orig_chord  = 3*(1-r_blade/R) + 1;

opt_twist   = optBlade(2,:);
opt_chord   = optBlade(1,:);



%% Plot

figure(101)
yyaxis left
plot(r_blade(1:end-1)/R,orig_cl,'-b','DisplayName','C_{L,original}')
hold on
grid on
plot(r_blade(1:end-1)/R,opt_cl,'--b','DisplayName','C_{L,opt}')

yyaxis right
plot(r_blade(1:end-1)/R,orig_cd,'-r','DisplayName','C_{D,original}')
plot(r_blade(1:end-1)/R,opt_cd,'--r','DisplayName','C_{D,opt}')
legend('show','Location','Northeast')


figure(102)
yyaxis left
plot(r_blade(1:end-1)/R,orig_cl,'-b','DisplayName','original C_L')
hold on
grid on
plot(r_blade(1:end-1)/R,opt_cl,'--b','DisplayName','optimal C_L')
ylabel('C_L')

yyaxis right
plot(r_blade(1:end-1)/R,orig_chord(1:end-1),'r-','DisplayName','original Chord')
plot(r_blade(1:end-1)/R,opt_chord(1:end-1),'r--','DisplayName','optimal Chord')
ylabel('Chord Distribution [m]')
legend('show','Location','Northeast')
xlabel('Radial Position r/R')


figure(200)

subplot(4,1,1)
    plot(r_blade(1:end-1)/R,orig_chord(1:end-1),'k-','DisplayName','original')
    hold on
    grid on
    plot(r_blade(1:end-1)/R,opt_chord(1:end-1),'b--','DisplayName','optimal')
    ylabel('Chord Distribution [m]','Interpreter','latex')
    legend('show','Location','Northeast')

subplot(4,1,2)
    plot(r_blade(1:end-1)/R,orig_twist(1:end-1),'k-','DisplayName','original')
    hold on
    grid on
    plot(r_blade(1:end-1)/R,opt_twist(1:end-1),'b--','DisplayName','optimal')
    ylabel('Twist Distribution $\theta [\deg]$','Interpreter','latex')
    legend('show','Location','Northeast')

subplot(4,1,3)
    plot(r_blade(1:end-1)/R,aoa_be,'k-','DisplayName','original')
    hold on
    grid on
    plot(r_blade(1:end-1)/R,aoa_be2,'b--','DisplayName','optimal')
    yline(optAngles(1),'g:','LineWidth',1.5,'DisplayName','\alpha for maximum C_L/C_D')
    yline(optAngles(2),'r:','LineWidth',1.5,'DisplayName','\alpha_{design} for optimum blade')
    ylabel('Angle of Attack $\alpha [\deg]$','Interpreter','latex')
    legend('show','Location','Southeast')

subplot(4,1,4)
    plot(r_blade(1:end-1)/R,orig_cl,'-k','DisplayName','original')
    hold on
    grid on
    plot(r_blade(1:end-1)/R,opt_cl,'--b','DisplayName','optimal')
    ylabel('Lift Coefficient [-]','Interpreter','latex')
    legend('show','Location','Southeast')

xlabel('Radial Position r/R','Interpreter','latex')