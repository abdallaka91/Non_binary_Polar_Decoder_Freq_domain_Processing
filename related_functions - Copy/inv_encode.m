function u = inv_encode(x1)
N=length(x1);
n=log2(N);

m1=nan(N,n+1);
m1(:,n+1)=x1;
for l=n:-1:1
    for t=0:N/2-1
        a=2*t-mod(t,2^(l-1))+1;
        b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
        tmp1 = [m1(a,l+1) m1(b,l+1)];
        m1(a,l)=bitxor(tmp1(1), tmp1(2));
        m1(b,l)=tmp1(2);
    end
end
u=m1(:,1);