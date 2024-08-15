n2=   %n2=2N
prheatmap=zeros(100,400);
parfor i=1:100
    for j=1:400
        prheatmap(i,j)=minpr(n2,(i-1)/100,(j-1)/100)%gives the matrix for the density plot for min PR/n2
    end
end
h=heatmap(prheatmap,'Colormap',flip(hot));
grid off
h.NodeChildren(3).YDir='normal';
XLabels = 0.01:0.01:1;
CustomXLabels = string(XLabels);
CustomXLabels(mod(XLabels,0.10) ~= 0) = " ";
YLabels = 0.01:0.01:4;
CustomYLabels = string(YLabels);
CustomYLabels(mod(YLabels,0.40) ~= 0) = " ";
xlabel('j_2/j') 
ylabel('g/j') 
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;


nk=1:350;
j0=zeros(length(nk),1);
j25=zeros(length(nk),1);
j5=zeros(length(nk),1);
j75=zeros(length(nk),1);
j1=zeros(length(nk),1);
parfor i=1:length(nk)
    j0(i)=minpr(2*nk(i),1.5,0);
    j25(i)=minpr(2*nk(i),1.5,0.25);
    j5(i)=minpr(2*nk(i),1.5,0.5);
    j75(i)=minpr(2*nk(i),1.5,0.75);
    j1(i)=minpr(2*nk(i),1.5,1);
end
hold on 
plot(nk,j0)
plot(nk,j0,'g')
plot(nk,j0,'y')
plot(nk,j0,'m')
plot(nk,j0,'r')
hold off
legend('j_2/j=0','j_2/j=0.25','j_2/j=0.5','j_2/j=0.75','j_2/j=1')
xlabel('2N') 
ylabel('min(PR/2N)')


n1=   %n1=2N
mgheatmap=zeros(100,400);
parfor i=1:100
    for j=1:400
        mgheatmap(i,j)=meangap(n1,(i-1)/100,(j-1)/100)%gives the matrix for the density plot for mean gap ratio
    end
end
h=heatmap(mgheatmap,'Colormap',flip(hot));
grid off
h.NodeChildren(3).YDir='normal';
XLabels = 0.01:0.01:1;
CustomXLabels = string(XLabels);
CustomXLabels(mod(XLabels,0.10) ~= 0) = " ";
YLabels = 0.01:0.01:4;
CustomYLabels = string(YLabels);
CustomYLabels(mod(YLabels,0.40) ~= 0) = " ";
xlabel('j_2/j') 
ylabel('g/j')
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomYLabels;

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

function prn = minpr(n,g,j2) %gives min(PR)/n for given parameters n, g, j2
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
