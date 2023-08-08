function res = calcM(L,KK)
    res=zeros(1,KK);
    for m=1:KK
        u=@(x) sin(m*pi/L*x).^2;
        res(m)=integral(u,0,L);
    end
end

