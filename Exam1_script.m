clear all
close all
clf

%% ----- Parameters-----%%

param.u=0.042*24;       %[m/day] - Sinking velocity of phytoplankton
param.D=10;             %[m^2/day] - Diffusion rate (for both plankton and nutrients)/ turbulent diffusion  
param.kphi=6e-10;       %[m^2/cell] - Specific light attenuation of phytoplankton
param.I0=600;           %[mmol photons/m^2 s^-1] - Incident light intensity
param.mu_max=0.94;      %[1/day]- Maximum specific growth rate of phytoplankton 
param.H_I=40;           %[µmol photons(I)/m2 s] - Half saturation constant of light limited growth
param.H_N= 0.0425;      %[mmol N/m^3] - Half saturation constant of nutrient limited growth
param.L=0.01*24;        %[1/day] - Specific loss rate in hours for phytoplankton
param.alpha=1*10^-9;    %[mmol N/(m^2s)] - Nutrient content of phytoplankton 
param.eps=0.005;        %[dimensionless] - Nutrient recycling coefficient
param.N_bottom = 50;    %[mmol N/m^3] - Bottom concentration of nutrients
param.w=5               %[m/day] - Sinking velocity of detritus
param.tau=0.1           %[1/day] - Detritus remineralisation rate

%System: 
param.n=50;             % Number of grid cells/points
param.Kbg=0.045;        %[1/m] - Background turbidity
depth=100;              %[m] - Depth of water column
param.dz=depth/param.n; %Width of section 
param.z=0.5*param.dz:param.dz:(depth-0.5*param.dz); %[1/m] - Grid definition



%% ----- Initital conditions ----- %% 

P0 = 2e6*exp(-(param.z-depth/4).^2/1000); 
N0 = param.N_bottom*exp(-(param.z-depth/1.8).^2/500); %normal distribution
D0=zeros(1,param.n);

PND=[P0,N0,D0]; %new vector with all initial conditions. 

%% ----- Run ODE45 -----%
[t,y] = ode45(@Exam1_function,[0:3100], PND, [], param);

%----- Labelling P, N, D and I -----%

P = y(:,1:param.n);
N = y(:,param.n+1:2*param.n);
D = y (:,2*param.n+1:end);
I = Exam1_calclight_function(P(end,:),t,param); 


%% ----- Make figures -----%

%Figure of 
figure(1)
subplot(4,1,1)
contourf(t,-param.z,P')
title('Phytoplankton')
ylabel('Depth(m)')
xlabel('Time(hours)')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
c = colorbar;
c.Label.String = 'Concentration of P [mmol/m^3]'
subplot(4,1,2)
title('Nutrients')
contourf(t,-param.z,N')
ylabel('Depth(m)')
xlabel('Time(hours)')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
c = colorbar;
c.Label.String = 'Concentration of N [mmol N/m^3]'
title('Nutrients')
subplot(4,1,3)
contourf(t,-param.z,D')
ylabel('Depth(m)')
xlabel('Time(hours)')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
c = colorbar;
c.Label.String = 'Concentration of D [mmol N/m^3]'
title('Detritus')
subplot(4,1,4)
contourf(t,-param.z,I)
ylabel('Depth[m]')
xlabel('Time[days]')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
c = colorbar;
c.Label.String = '[mmol photons/m^2 s^-1]'
title('Light')

%Convergence(t3100)WINTER:
figure(2)
subplot(1,3,1)
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
plot(P(end,1:end),-param.z,'Linewidth',2)
title('Phytoplankton, winter')
ylabel('Depth[m]')
xlabel('Conc. of P [mmol N/m^3]')
subplot(1,3,2)
plot(N(end,1:end),-param.z,'Linewidth',2)
title('Nutrients, winter')
ylabel('Depth[m]')
xlabel('Conc. of N [mmol N/m^3]')
subplot(1,3,3)
plot(D(end,1:end),-param.z,'Linewidth',2)
title('Detritus, winter')
ylabel('Depth[m]')
xlabel('Conc. of D [mmol N/m^3]')

%Convergence(t2918)SUMMER:
figure(3)
subplot(1,3,1)
plot(P(2918,1:end),-param.z,'Linewidth',2)
title('Phytoplankton, summer')
ylabel('Depth[m]')
xlabel('Conc. of P [mmol N/m^3]')
subplot(1,3,2)
plot(N(2918,1:end),-param.z,'Linewidth',2)
title('Nutrients, summer')
ylabel('Depth[m]')
xlabel('Conc. of N [mmol N/m^3]')
subplot(1,3,3)
plot(D(2918,1:end),-param.z,'Linewidth',2)
title('Detritus, summer')
ylabel('Depth[m]')
xlabel('Conc. of D [mmol N/m^3]')


%Limiting factors of P- WINTER
figure(4)
subplot(1,3,1)
plot(P(end,1:end),-param.z,'Linewidth',2)
title('Phytoplankton')
ylabel('Depth[m]')
xlabel('Conc. of P [mmol N/m^3]')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
yline(-18)
subplot(1,3,2)
plot(N(end,1:end),-param.z,'Linewidth',2)
title('Nutrients')
ylabel('Depth[m]')
xlabel('Conc. of N [mmol N/m^3]')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
yline(-18)
subplot(1,3,3)
plot(I(1:end,end),-param.z,'Color','#0072BD','Linewidth',2)
title('Light')
ylabel('Depth[m]')
xlabel('I [mmol photons/m^2 s^-1]')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
yline(-18)

% D with different w. Remember to change w
figure(5)
plot(D(2918,1:end),-param.z,'Linewidth',2)
title('Detritus, w=x')
ylabel('Depth[m]')
xlabel('Conc. of D [mmol N/m^3]')

%figure of different eps. Remeber to change!!
figure(6)
subplot(3,1,1)
contourf(t,-param.z,P')
title('Phytoplankton')
ylabel('Depth(m)')
xlabel('Time(hours)')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
c = colorbar;
c.Label.String = 'Concentration of P [mmol/m^3]'
subplot(3,1,2)
title('Nutrients, eps=') %remeber to change!
contourf(t,-param.z,N')
ylabel('Depth(m)')
xlabel('Time(hours)')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
c = colorbar;
c.Label.String = 'Concentration of N [mmol N/m^3]'
title('Nutrients')
subplot(3,1,3)
contourf(t,-param.z,D')
ylabel('Depth(m)')
xlabel('Time(hours)')
set(gca,'FontName','Times New Roman','FontSize',13,'Ycolor','k')
c = colorbar;
c.Label.String = 'Concentration of D [mmol N/m^3]'
title('Detritus')

