function dydt = Exam1_function(t,y,param)  

P=y(1:param.n);
N=y(param.n+1:2*param.n);
D=y(2*param.n+1:end);

%% -----Diffussion-Advection for Plankton-----%%

%surface boundary fluxes:
Ja_P(1)=0; %no advection at the surface 
Jd_P(1)=0; %no difussion at the surface
%bottom boundary fluxes:
Ja_P(param.n+1)=0; %no advection at the bottom
Jd_P(param.n+1)=0; %no difussion at the bottom

for i=2:param.n 
    Ja_P(i)=param.u*P(i-1); %Advection 
    Jd_P(i)=-param.D*(P(i)-P(i-1))/param.dz; % Diffusion
end

J_plankton=Ja_P+Jd_P; %added together to get the whole water column. 

for i=1:param.n 
    dPdt(i)=-(J_plankton(i+1)-J_plankton(i))/param.dz; 
end
%% -----Diffussion-Advection for Nutrients-----%%

for i=2:param.n
    Ja_N(i)=0; %nutrients don't have any sinking velocity therefore advection=0
    Jd_N(i)=-param.D.*(N(i)-N(i-1))/param.dz;
end

%surface boundary fluxes:
Ja_N(1)=0;
Jd_N(1)=0;
%bottom boundary fluxes
Ja_N(param.n+1)=0; 
Jd_N(param.n+1)=-param.D*((param.N_bottom-N(param.n))/param.dz); 

J_nutrients=Ja_N+Jd_N; %added together to get the whole water column. 

for i=1:param.n
    dNdt(i)=-(J_nutrients(i+1)-J_nutrients(i))/param.dz;
end


%% -----Diffussion-Advection for Detrius-----%%

%surface boundary fluxes:
Ja_D(1)=0; 
Jd_D(1)=0;
%bottom boundary fluxes
Ja_D(param.n+1)=0; 
Jd_D(param.n+1)=0; 

for i=2:param.n
    Ja_D(i)=param.w*D(i-1); 
    Jd_D(i)=-param.D*(D(i)-D(i-1))/param.dz;
end

J_detritus = Ja_D + Jd_D; %added together to get the whole water column. 

for i=1:param.n
    dDdt(i)=-(J_detritus(i+1)-J_detritus(i))/param.dz;
end

%% ----- LIGHT -----%%

I=Exam1_calclight_function(P',t,param); %Kalder den hen til lysfuntionen "calclight_function"


%% ----- GROWTH -----%%
mu=param.mu_max.*min(I./(param.H_I+I),N./(param.H_N+N)); %min chooses the minimum elements of the arrays. Used for Liebig's law.

%% ----Update derivaties---%% 

dPdt=dPdt'+mu.*P-param.L.*P;
dNdt=dNdt'-param.alpha*mu.*P+param.tau*D+param.eps*param.alpha*param.L.*P;
dDdt=dDdt'-param.tau*D+param.alpha.*param.L.*P;

dydt=[dPdt;dNdt;dDdt];
end