clear all
close all
clf
%%
dz = [0.2 0.8 1.4 2.0 2.6 3.2 6.2];

for i = 1:length(dz)
    [t,z,P] = Gridsensitivity_function(dz(i));
    plot(P(end,:), -z,'--.','Linewidth',1)
    drawnow
    hold on
end
%%
legend('dz=0.2','dz=0.8','dz=1.4','dz=2','dz=2.6','dz=3.2','dz=6.2')
xlabel('Concentration of Phytoplankton [mmol/m^3]')
ylabel('Depth (m)')
