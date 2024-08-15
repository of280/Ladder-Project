j2= ; %enter value for j2/j
qd=-pi:0.01:pi;
evcont=zeros(length(qd),2);
devcont=zeros(length(qd),2);
for i=1:length(qd)
    evcont(i,1)=evup(1.5,j2,qd(i));
    evcont(i,2)=evdn(1.5,j2,qd(i));
    devcont(i,1)=dvup(1.5,j2,qd(i));
    devcont(i,2)=dvdn(1.5,j2,qd(i));
end

hold on%plots the group velocities for your value fo j_2/j at g=1.5j
plot(qd,devcont(:,1),'b')
plot(qd,devcont(:,2),'m')
yline(0)
xlim([-pi,pi])
xlabel('qd') 
ylabel('v_{q,\tau}/jd')
hold off

hold on%plots the eigenvalues for your value of j_2/j at g=1.5j
plot(qd,evcont(:,1),'b')
plot(qd,evcont(:,2),'m')
xlabel('qd') 
ylabel('(\omega-\omega_0)/j'
xlim([-pi,pi]))
hold off
%

nk= ; %nk=2N
prs=minpr(nk,1.5,j2); %calculates eigenvalues and associated PR/2N for your j2/j and nk=2N
[~,idx]=min(prs(:,2)); %index of the minPR eigenvalue
hold on
scatter(prs(:,2),prs(:,1),'.')
yline([max(evcont(:,1)),min(evcont(:,1))],'b')
yline([max(evcont(:,2)),min(evcont(:,2))],'m')
yline(prs(idx,1),'black')
xlabel('PR/2N') 
ylabel('(\omega-\omega_0)/j')
hold off


function rr = evdn(g,j2,x)%lower eigenvalue
    rr=2*cos(x)+j2*cos(2*x)-sqrt(2)*sqrt(2*g^2+j2^2+(j2^2)*cos(4*x))/2;
end
function rr = evup(g,j2,x)%upper eigenvalue
    rr=2*cos(x)+j2*cos(2*x)+sqrt(2)*sqrt(2*g^2+j2^2+(j2^2)*cos(4*x))/2;
end
function rr = dvup(g,j2,x)%upper eigenvalue first derivative
    rr=-2*sin(x)-2*j2*sin(2*x)-sqrt(2)*(j2^2)*sin(4*x)/(2*g^2+j2^2+(j2^2)*cos(4*x))^(1/2);
end
function rr = dvdn(g,j2,x)%lower eigenvalue derivative
    rr=-2*sin(x)-2*j2*sin(2*x)+sqrt(2)*(j2^2)*sin(4*x)/(2*g^2+j2^2+(j2^2)*cos(4*x))^(1/2);
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
    [v,d]=eig(A);
    ipr=zeros(n,1);
    ev=zeros(n,1);
    for i=1:n
         ipr(i)=sum((v(:,i).^4)).^(-1);
         ev(i)=d(i,i);
    end
    prn=[ev,ipr/n];

end
