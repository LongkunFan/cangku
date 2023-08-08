function res=NumSol(L,W,V,R,R1,R2,D,T,B,tau_q,tau_p)
    k=0.01;
    h=0.05;
    L=floor(L/h)*h;
    W=floor(W/h)*h;
    T=floor(T/k)*k;
    xx=0:h:L;
    yy=0:h:W;
    tt=0:k:T;
    N=T/k+1;
    M=floor(L/h)+1;
    r=h/k;
    res=zeros(M,M,N+1);%N=1仅为方便处理第一层纽曼边界条件,从N=2开始才是结果
    CC=tau_q*r^2;
    AA=4+r*h+CC;
    BB=r*h+2*CC;
    myI=-eye(M-2);
    myB=eye(M-2)*AA;
    A=zeros((M-2)^2,(M-2)^2);
    for i=1:(M-2)
        A((i-1)*(M-2)+1:i*(M-2),(i-1)*(M-2)+1:i*(M-2))=myB;
    end
    for i=1:M-3
        A((i-1)*(M-2)+1:i*(M-2),i*(M-2)+1:(i+1)*(M-2))=myI;
        A(i*(M-2)+1:(i+1)*(M-2),(i-1)*(M-2)+1:i*(M-2))=myI;
    end
    K=zeros((M-2)^2,1);
    
    for jj=2:N%实际计算时刻
        j=jj+1;%当前计算结果res存储层数
        K(1)=res(1,2,j)+res(2,1,j);
        K((M-3)*(M-2)+1)=res(1,M-1,j)+res(2,M,j);%(dirichlet_func(a,a+h)+dirichlet_func(a+h,a));
        K(M-2)=res(M-1,1,j)+res(M,2,j);%(dirichlet_func(b-h,b)+dirichlet_func(b,b-h));
        K((M-3)*(M-2)+M-2)=res(M-1,M,j)+res(M,M-1,j);%(dirichlet_func(b-h,a)+dirichlet_func(b,a+h));
        for i=2:M-3
            K(i)=res(i+1,1,j);%dirichlet_func(a+i*h,b);
            K((M-3)*(M-2)+i)=res(i+1,M,j);%dirichlet_func(a+i*h,a);
            K((i-1)*(M-2)+M-2)=res(M,i+1,j);%dirichlet_func(b,b-i*h);
            K((i-1)*(M-2)+1)=res(1,i+1,j);%dirichlet_func(a,b-i*h);
        end
        for m=1:M-2
            for n=1:M-2
                K((m-1)*(M-2)+n)=K((m-1)*(M-2)+n)+...
                    ((1-B)*(res(m,n+1,jj)+res(m+1,n,jj)+res(m+2,n+1,jj)+res(m+1,n+2,jj)-4*res(m+1,n+1,jj))+...
                    BB*res(m+1,n+1,jj)-CC*res(m+1,n+1,jj-1)+myQ(D,xx(m+1),yy(n+1),tt(jj),tau_p,R,V,W));
            end
        end
        res(2:M-1,2:M-1,j)=reshape(A\K,M-2,M-2)';
    end
    
end