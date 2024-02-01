function prn = minipr(n,g,j2) %gives min(IPR)/n for given parameters n, g, j2
    e = zeros(n,1);
    f = zeros(n,1);
    h = zeros(n,1);
    for i = 1:2:n-1
        e(i)=g;
        e(i+1)=0;
        f(i)=1;
        f(i+1)=1;
        h(i)=j2;
        h(i+1)=0;
    end
    A = full(spdiags([h f e flip(e) f h],[-4 -2 -1 1 2 4],n,n));
    [v,~]=eig(A);
    ipr=zeros(n,1);
    for i=1:n
         ipr(i)=sum((v(:,i).^4)).^(-1);
    end
    prn=min(ipr)/n;
end