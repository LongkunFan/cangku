function res = calcN(W,KK)
    res=zeros(1,KK);
    for n=1:KK
        u=@(y) sin(n*pi*y/W).^2;
        res(n)=integral(u,0,W);
    end
end

