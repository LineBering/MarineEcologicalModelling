function I=Exam1_calclight_function(P,t,param)       
int=cumsum(P.*param.dz)*param.kphi; %cumsum is the cumulative sum of elements.
I=(param.I0*(1+cos((2*pi*t)/365)))*exp(-param.Kbg*param.z-int); 
I=I';
end