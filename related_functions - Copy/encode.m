function x1 = encode(u)
N=length(u);
n=log2(N);

m1=nan(N,n+1);
m1(:,1)=u;
for l=1:n
    for t=0:N/2-1
        a=2*t-mod(t,2^(l-1))+1;
        b=2^(l-1)+2*t-mod(t,2^(l-1))+1;
        tmp1 = [m1(a,l-1+1) m1(b,l-1+1)];
        m1(a,l-1+2)=bitxor(tmp1(1), tmp1(2));
        m1(b,l-1+2)=tmp1(2);
    end
end
x1=m1(:,end);