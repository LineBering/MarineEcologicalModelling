
function [t,z,P] = Gridsensitivity_function(dz) % dz is the grid spacing

% Parameters:
param.D = 10; 
param.u = 0.042*24; 

param.depth = 100; 
param.dz = dz; 
param.z = param.dz/2:param.dz:(param.depth-param.dz/2);
param.n = length(param.z); 
P0 = 2e6*exp(-(param.z-param.depth/4).^2/1000);

[t, P] = ode45(@adv_dif_P, [0:200], P0, [], param);
z = param.z;

    function dPdt = adv_dif_P(t,P,param)
        
        ix = 2:(param.n);
        Jadv(ix) = param.u*P(ix-1);
        Jadv(1) = 0;
        Jadv(param.n+1) = 0; 
     
        Jdiff(ix) = -param.D*(P(ix)-P(ix-1))/param.dz;
        Jdiff(1) = 0;
        Jdiff(param.n+1) = 0; 
        J = Jadv + Jdiff;
        dPdt = -(J(2:(param.n+1))-J(1:param.n))/param.dz;

        dPdt = dPdt';
    end

end


