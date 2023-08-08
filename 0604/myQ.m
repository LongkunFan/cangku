function res=myQ(D,X,Y,tau,tau_p,R,V,W)
res=D*exp(-(tau/tau_p).^2-2*((X-V*tau).^2+(Y-W/2).^2)/R^2)*(1+2*V*(X-V*tau)/R^2-tau/tau_p^2);
end