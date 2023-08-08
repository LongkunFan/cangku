function res = calcFmn(L,V,R,T,tau_p,KK,mybeta1,mybeta2)
syms t;    
res=zeros(KK,KK);
    for m=1:KK
        for n=1:KK

            u=@(x,t) sin(m*pi/L*x).*(exp(-2*((x-V*t)/R).^2+(mybeta2(m,n)-mybeta1(m,n))*(T-t))...
                -exp(-2*((x-V*t)/R).^2-(mybeta1(m,n)+mybeta2(m,n))*(T-t)))...
                .*(1+(2*V*(x-V*t)/R^2))/2;
            
% %             u=@(x,t) sin(m*pi/L*x).*(exp(-(t/tau_p).^2-2*((x-V*t)/R).^2+(mybeta2(m,n)-mybeta1(m,n))*(T-t))...
% %                 -exp(-(t/tau_p).^2-2*((x-V*t)/R).^2-(mybeta1(m,n)+mybeta2(m,n))*(T-t)))...
% %                 .*(1-t/tau_p^2+(2*V*(x-V*t)/R^2))/2;

% %             u=@(x,t) sin(m*pi/L*x).*(exp(-2*((x-V*t)/R).^2+(mybeta2(m,n)-mybeta1(m,n))*(T-t))...
% %                 -exp(-2*((x-V*t)/R).^2-(mybeta1(m,n)+mybeta2(m,n))*(T-t)))...
% %                 .*(1+(2*V*(x-V*t)/R^2))/2;
            f=integral2(u,0,L,0,T);
            res(m,n)=f;
        end
    end
end

