function val=calc(X,Y,L,W,mybeta2,myMm,myNn,myFn,myFmn,mycccs,KK)
%T «º∆À„ ±øÃ
val=0;
for m=1:KK
    for n=1:KK
        Am=m*pi/L;Bn=n*pi/W;
        %t=sin(Am*X)*sin(Bn*Y)*mycccs(m,n);
        t=sin(Am*X)*sin(Bn*Y);
        val=val+t*(myFn(n)*myFmn(m,n))/(myMm*myNn*mybeta2(m,n));
    end
end
%{
mybeta2=conj(mybeta2);
for m=1:KK
    for n=1:KK
        Am=m*pi/L;Bn=n*pi/W;
        t=sin(Am*X)*sin(Bn*Y)*mycccs(m,n);
        val=val+t*(myFn(n)*myFmn(m,n))/(myMm*myNn*mybeta2(m,n));
    end
end
val=val/2;
%}

end

