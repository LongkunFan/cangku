function res = calcLambda2(L,W,KK)
    res=zeros(KK,KK);
    for i=1:KK
        res(i,:)=(i*pi/L)^2+((1:KK)*pi/W).^2;
    end
end

