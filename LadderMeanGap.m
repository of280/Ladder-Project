function mg = meangap(n,g,j2)  %gives the mean gap ratio for given values of n, g, j2
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
    ev=eig(A);
    spc=zeros(length(ev)-2,1);
    gr=zeros(length(ev)-2,1);
    for i=1:length(spc)
        spc(i)=(ev(i+2)-ev(i+1))/(ev(i+1)-ev(i));
        gr(i)=min(spc(i),1/spc(i));
    end
    mg=mean(gr);
end