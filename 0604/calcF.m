function res = calcF(W,R,KK)
%{
    res=zeros(1,KK);
    for n=1:KK
        u=@(y) sin(n*pi*y/W).*exp(-2*(y-W/2).^2/R^2);
        res(n)=integral(u,0,W);
    end
  %}
N=128;
f=@(x) exp(-2*(x-W/2).^2/R^2);
U=zeros(1,N);
for i=2:N/2
    U(i)=f(2*W/(N)*(i-1));
    U(N-i+2)=-f(2*W/(N)*(i-1));
end
U(1)=f(0);U(N/2+1)=-f(W);
V=fft(U);
FF2=V(2:N/2+1);
FF2(N/2)=FF2(N/2)/2;
FF2=-imag(FF2*2/N*(W/2));
res=FF2(1:KK);
end

