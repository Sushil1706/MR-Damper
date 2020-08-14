i = 1;
r = 30/1000;
rc = 10/1000;
lp = 20/1000;
t = 8/1000;
l = 50/1000;
wc = 20/1000;
n=0.112;
y1=0.4;
nc=150;
u0=4*(pi)*10^(-7);
um=6;
u=600;
c0=0.3;
c1=0.42;
c2=-0.00116;
c3=1.0513*10^-6;


re=85.1*10^(-3);
lw=12;
M=[];
N=[];
V=[];
K=[];
J=[];
H=[];

for g = 0.001:0.0001:0.005
    rd=r+g/2;
    l1 = r-rc/2;
    l7 = l1;
    l2 =  g;
    l6 = l2;
    l3 = t/2;
    l5 =l3;
    l4 = l-lp;
    l8 = l4;
    a1 = 2*pi*lp*(r-rc/4);
    a7 = a1;
    a2 = 2*pi*lp*(r+g);
    a6 = a2;
    a3 = 2*pi*lp*(r+g+t/4);
    a5 = a3;
    a4 = pi*((r+g+t)^2 -(r+g+t/2)^2);
    a8 = pi*rc^2;
    hmr =((nc*i)/((2*g)+((2*l1*um*a2)/(u*a1))+((2*l3*um*a2)/(u*a3)) + ((l4*um*a2)/(u*a4))+ ((l8*um*a2)/(u*a8))));
    bmr = u0*um*hmr;
    ap=pi*(r+g)^2;
    as=pi*(r)^2;
    flux = a2*bmr;
    lin=nc*(flux/i);
    rw = lw*re;
    if y1==0
         sgn(y1) = 0;
    elseif y1>0
         sgn('y1') =1;
    elseif y1<0
         sgn(y1)= -1;
    end
    
    ty = c0+c1*hmr+c2*hmr^2+c3*hmr^3;
    q=y1*(ap-as);
    c=2.07+12*q*n/(12*q*n+0.8*pi*rd*g^2*ty);

    cv=6*n*l*(ap-as)^2/(pi*rd*g^3);
    fm=(ap-as)*2*c*lp*ty/g;
    fv=cv*y1;
    fd=fv+(fm*sgn('y1'));
    Tk=lin/rw;
    M=[M fd(1)];
    dy=(fv+fm)/fv;
    N=[N dy];
    V=[V g];
    K=[K ty];
    J=[J Tk];
    H =[H hmr];
end
figure(1)
plot(V,N)

% figure(2)
% plot(V,M)
% 
% figure(3)
% plot(V,K)
% 
% figure(4)
% plot(V,J)