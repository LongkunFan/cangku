function res = cccs(A,B,R1,R2)
K=1:35;
if A~=B
    res=sum((-1).^K.*...
        (A*R1.*besselj(2*K,B*R1).*besselj(2*K-1,A*R1)...
        -B*R1.*besselj(2*K,A*R1).*besselj(2*K-1,B*R1)...
        -A*R2.*besselj(2*K,B*R2).*besselj(2*K-1,A*R2)...
        +B*R2.*besselj(2*K,A*R2).*besselj(2*K-1,B*R2)));
    res=(4*pi/(A^2-B^2))*(res+0.5*(...
        -A*R1.*besselj(0,B*R1).*besselj(1,A*R1)...
        +B*R1.*besselj(0,A*R1).*besselj(1,B*R1)...
        +A*R2.*besselj(0,B*R2).*besselj(1,A*R2)...
        -B*R2.*besselj(0,A*R2).*besselj(1,B*R2)));
else
    res=(1/(2*A))*sum((-1).^K.*(-A*R1^2*besselj(2*K,A*R1).^2+A*R2^2*besselj(2*K,A*R2).^2 ...
    +4*K.*R1.*besselj(2*K,A*R1).*besselj(2*K-1,A*R1)...
    -A*R1^2*besselj(2*K-1,A*R1).^2-4*K.*R2.*besselj(2*K,A*R2).*besselj(2*K-1,A*R2)...
    +A*R2^2*besselj(2*K-1,A*R2).^2));
    res=res+(-R1^2*(besselj(0,A*R1)^2+besselj(1,A*R1)^2)+R2^2*(besselj(0,A*R2)^2+besselj(1,A*R2)^2))/2*0.5;
    res=4*pi*res;
end
end

