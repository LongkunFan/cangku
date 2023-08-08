function [res,mybeta1,mybeta2,myMm,myNn,myFn,myFmn,mycccs]=main(L,W,V,R,R1,R2,D,T,B,tau_p)
N=300;%x
M=300;%y
res=zeros(N+1,M+1);
%KKÎª½Ø¶Ï
KK=40;
mylambda2=calcLambda2(L,W,KK);
mybeta1=1+(B/2)*mylambda2;
mybeta2=sqrt((2+B*mylambda2).^2-4*mylambda2)/2;
myMm=L/2;
myNn=W/2;
myFn=calcF(W,R,KK);
myFmn=calcFmn(L,V,R,T,tau_p,KK,mybeta1,mybeta2);
mycccs=zeros(KK,KK);
for m=1:KK
    for n=1:KK
        Am=m*pi/L;Bn=n*pi/W;
        mycccs(m,n)=cccs(Am,Bn,R1,R2);
    end
end

for i=1:N+1
    X=(i-1)*(L/N);
    for j=1:M+1
        Y=(j-1)*(W/M);
        res(i,j)=D*calc(X,Y,L,W,mybeta2,myMm,myNn,myFn,myFmn,mycccs,KK);
    end
end
end