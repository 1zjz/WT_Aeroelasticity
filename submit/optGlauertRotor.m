clear;clc;close all;
%% Glauert Optimum Rotor
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
UsePrandtl = 1;       % Flag for using Prandtl correction
% Calculate radial distribution (linearly spaced)
r_blade = linspace(r0, R, N);
% Calculate twist distribution
twist_blade = 14*(1-r_blade/R);
% Calculate chord distribution 
chord_blade = 3*(1-r_blade/R) + 1;

% Parameters task 1
U0 = 10;                     % Downstream wind velocity [m/s]
TSR = 8;             % [6 8 10]; % range of TSR
file_name = 'DU95W180.csv'; % cl,cd,cm,alpha data for airfoil
omega = U0*TSR/R;           % Rotor speeds for TSR
CTlim = 0.75;               % CT limit parameter

%% Load Polar Curves of the aerofoil and plot them
polar_data = table2array(readtable(file_name)); 
AOA = polar_data(:,1); %[deg]
CL  = polar_data(:,2);
CD  = polar_data(:,3);
CM  = polar_data(:,4);

% % % % % %% Design AoA
% % % % % cl_cd               = CL ./ CD;
% % % % % [tempmax,tmaxidx]   = max(cl_cd);
% % % % % opt_aoa            = AOA(tmaxidx);
% % % % % opt_cl             = CL(tmaxidx);
% % % % % opt_cd             = CD(tmaxidx);
% % % % % 
% % % % % [sort_clcd,sort_idx] = sort(cl_cd,'descend');
% % % % % sort_AOA = AOA(sort_idx);
% % % % % 
% % % % % % Preallocation
% % % % % theta_glauert = zeros(1,length(r_blade)); % glauert optimum pitch
% % % % % c_glauert     = zeros(1,length(r_blade)); % glauert optimum chord
% % % % % CT_glauert    = zeros(1,length(AOA)); % glauert CT // check >= 0.75
% % % % % cn            = zeros(1,length(r_blade)); % normal force coefficient
% % % % % 
% % % % % for ii = 1:length(AOA)
% % % % %     
% % % % %     alpha = AOA(ii);              % extract aoa
% % % % %     cl  = interp1(AOA,CL,alpha);        % 
% % % % %     cd  = interp1(AOA,CD,alpha);
% % % % %     
% % % % %     for jj = 1:N
% % % % %         
% % % % %         
% % % % %         
% % % % %         x = TSR*r_blade(jj)/R;
% % % % %         
% % % % %         % Compute a from third order polynom \cite WTTA BEM2 class @DTU
% % % % %         syms aa
% % % % %         eqn = 16*aa^3 - 24*aa^2 + aa*(9-3*x^2)-1+x^2 == 0;
% % % % %         aa_sol = vpasolve(eqn,aa);
% % % % %         atemp = aa_sol(aa_sol > 0);
% % % % %         atemp = atemp(atemp < 1);
% % % % %         a = double(atemp); clear atemp;
% % % % %         
% % % % %         ap = (1-3*a)/(4*a-1);
% % % % %         phi = rad2deg(atan((1-a)/((1+ap)*x))); % Inflow angle
% % % % %         theta_glauert(jj) = phi-alpha; % resulting Pitch
% % % % %         
% % % % %         % Prandtl Top and Root Loss Correction
% % % % %         %F = 2/pi * acos((exp(-B/2*(R-r_blade(jj))/(r_blade(jj)*sind(phi))))); % Tip loss factor
% % % % %         [F,~,~] = PrandtlTipRootCorrection(r_blade/R, 0/R,r_blade(end)/R,TSR,B,a);
% % % % %         %[F,~,~] = PrandtlTipRootCorrection(r_blade, 0,r_blade(end),TSR,B,a);
% % % % %         %[F,~,~] = PrandtlTipRootCorrection(r_blade/R, 0.2,r_blade(end)/R,TSR,B,a);
% % % % %         %[F,~,~] = PrandtlTipRootCorrection(r_blade, r_blade(1),r_blade(end),TSR,B,a);
% % % % %         cn(jj) = cl*cosd(phi) + cd*sind(phi); % Force c in normal for CT
% % % % % 
% % % % %         c_glauert(jj) = R*8*pi*F(jj)*a*x*sind(phi)^2/((1-a)*B*TSR*cn(jj));
% % % % %     end
% % % % %     
% % % % %     % Flag when the target of CT = 0.75 is met
% % % % %     CT_glauert(ii) = trapz(r_blade/R,cn(:));
% % % % % end
% % % % % % Save Values as .mat
% % % % % data_optGlauert = [AOA,cl_cd,CT_glauert'];
% % % % % save('optGlauert.mat','data_optGlauert');


%% Load optGlauert.mat
% Skip commented for-loop by loading this file
load('optGlauert.mat');
AOA         = data_optGlauert(:,1);
cl_cd       = data_optGlauert(:,2);
CT_glauert  = data_optGlauert(:,3);

%% Plots
figure(101)
yyaxis left
    plot(AOA,cl_cd,'b-','DisplayName','C_L / C_D')
        ylabel('$C_{L} / C_{D}$ [-]', 'Interpreter','latex')
yyaxis right
    plot(AOA,CT_glauert,'r--','DisplayName','C_T')
    yline(CTlim,'r:','LineWidth',1.75,'DisplayName','C_{T,target}')
        ylabel('$C_{T}$', 'Interpreter','latex')
xlabel('Angle of Attack [deg]')
legend('show','Location','Southeast')

%% Determine Design Angle of Attack
idx_low     = find(floor(AOA) == 0);
idx_top     = find(floor(AOA) == 20);
designrange = idx_low:idx_top;

aoa_design = interp1(CT_glauert(designrange),AOA(designrange),CTlim);

%% Glauert Optimum Rotor

cl  = interp1(AOA,CL,aoa_design);        % 
cd  = interp1(AOA,CD,aoa_design);

for jj = 1:N

    x = TSR*r_blade(jj)/R;

    % Compute a from third order polynom \cite WTTA BEM2 class @DTU
    syms aa
    eqn = 16*aa^3 - 24*aa^2 + aa*(9-3*x^2)-1+x^2 == 0;
    aa_sol = vpasolve(eqn,aa);
    atemp = aa_sol(aa_sol > 0);
    atemp = atemp(atemp < 1);
    a = double(atemp); clear atemp;

    ap = (1-3*a)/(4*a-1);
    phi = rad2deg(atan((1-a)/((1+ap)*x))); % Inflow angle
    theta_glauert(jj) = phi - aoa_design; % resulting Pitch

    % Prandtl Top and Root Loss Correction
    %F = 2/pi * acos((exp(-B/2*(R-r_blade(jj))/(r_blade(jj)*sind(phi))))); % Tip loss factor
    [F,~,~] = PrandtlTipRootCorrection(r_blade/R, 0/R,r_blade(end)/R,TSR,B,a);
    %[F,~,~] = PrandtlTipRootCorrection(r_blade, 0,r_blade(end),TSR,B,a);
    %[F,~,~] = PrandtlTipRootCorrection(r_blade/R, 0.2,r_blade(end)/R,TSR,B,a);
    %[F,~,~] = PrandtlTipRootCorrection(r_blade, r_blade(1),r_blade(end),TSR,B,a);
    cn(jj) = cl*cosd(phi) + cd*sind(phi); % Force c in normal for CT

    c_glauert(jj) = R*8*pi*F(jj)*a*x*sind(phi)^2/((1-a)*B*TSR*cn(jj));
end


%Substract the minimum twist to set = 0 at tip
twist_glauert = theta_glauert - theta_glauert(end);

%% Save optimum blade
optBlade = [c_glauert;twist_glauert];
save('optBlade.mat','optBlade');


%% Plot Twist and Chord Distribution

figure(102)
subplot(2,2,1:2)
plot(r_blade/R,twist_blade,'.k-','DisplayName','original Twist')
hold on
grid on
plot(r_blade/R,twist_glauert,'.b--','DisplayName','optimal Twist (Glauert Rotor)')
ylabel('Twist angle [deg]')
xlabel('r/R')
legend('show','Location','Northeast')

subplot(2,2,3:4)
plot(r_blade/R,chord_blade,'.k-','DisplayName','original Chord')
hold on
grid on
plot(r_blade/R,c_glauert,'.b--','DisplayName','optimal Chord (Glauert Rotor)')
ylabel('Chord [m]')
xlabel('r/R')
legend('show','Location','Northeast')

sgtitle(['Glauert Optimal Rotor at \alpha_{design} = ',num2str(round(aoa_design,4)), ' deg'])







%% Operational Point

clcd_opt    = interp1(AOA,cl_cd,aoa_design);
clcd_max    = max(cl_cd);
aoa_max     = interp1(cl_cd,AOA,clcd_max);


figure(103)
    plot(AOA,cl_cd,'k-','DisplayName','C_L / C_D')
    hold on
    grid on
    plot(aoa_design,clcd_opt,'r.','MarkerSize',15,'DisplayName','Operational Point: Design')
    plot(aoa_max,clcd_max,'g.','MarkerSize',15,'DisplayName','Operational Point: Max')
        ylabel('C_L / C_D')
        xlabel('Angle of Attack [deg]')
legend('show','Location','Southeast')

optAngles = [aoa_max,aoa_design];
save('optAngles.mat',"optAngles");




