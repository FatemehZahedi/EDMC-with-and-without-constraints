clear
clc
[n1,d1,n2,d2]=Inputsys(1);
Gs1 = tf(n1,d1);
Ts=0.1;
Gd1 = c2d(Gs1,Ts,'zoh');
[num1,den1]=tfdata(Gd1,'v');
Gs2 = tf(n2,d2);
Gd2 = c2d(Gs2,Ts,'zoh');
[num2,den2]=tfdata(Gd2,'v');
sys_info = stepinfo(Gd1);
ts1 = sys_info.SettlingTime;
tr1=sys_info.RiseTime; 
sys_info = stepinfo(Gd2);
ts2 = sys_info.SettlingTime;
tr2=sys_info.RiseTime;
t=1:Ts:30;
P1=floor(tr1/Ts);
P2=floor(tr2/Ts);
P=max(P1,P2);
M=P-2;
%-----------------------------------------------------------------------
gamma =1;
alpha=0.5;
betta=0.4;
delta=0.01;
% gain_DC=(num1(1)+num1(2)+num1(3))/(den1(1)+den1(2)+den1(3));
% gain_DC2=(num2(1)+num2(2)+num2(3))/(den2(1)+den2(2)+den2(3));
% Q = eye(P);
% R1 =((gain_DC2/gain_DC)^2)*gamma*gain_DC^2*eye(M);
% R2=gamma*gain_DC2^2*eye(M);
% R=[R1 zeros(M); zeros(M) R2];
%-----------------------------------------------------------------------
x01=0.0882;
x02=441.2;
%-----------------------------------------------------------------------
up1=0.4;
up2=-0.2;
low1=0.1;
low2=-0.7;
% up1=0.4;
% up2=0.2;
% low1=-0.5;
% low2=-0.9;
yre=[];
u_1=[];
u_2=[];
x_1=[];
yeli=[];
ym=[];
y=0;

Y_d=zeros(P,length(t));
Y_past=zeros(P,length(t));
Y_m=zeros(P,length(t));
Y_eli=zeros(P,length(t));
Y_nm=zeros(P,length(t));
Y_per=zeros(P,length(t));

U1_ = zeros(P,length(t));
U2_ = zeros(P,length(t));
U_=[U1_; U2_];

U1per = zeros(P,length(t));
U2per = zeros(P,length(t));
Uper=[U1per; U2per];

U1nm = zeros(P,length(t));
U2nm = zeros(P,length(t));

dU1=zeros(M,length(t));
dU2=zeros(M,length(t));
dU=[dU1;dU2];

dext=zeros(1,length(t));
Dext=zeros(P,length(t));
Dnl=zeros(P,length(t));
E=zeros(P,length(t));

g=zeros(P,length(t));

g1=zeros(P,length(t));
g2=zeros(P,length(t));
Y_per1=zeros(P,length(t));
Y_per2=zeros(P,length(t));
Y_past1=zeros(P,length(t));
Y_past2=zeros(P,length(t));

U1=dU1;
U2=dU2;

dU01=dU;
dU01(1,1)=0.001;
dU01(M+1,1)=0.001;

n=zeros(1,length(t));

noise=zeros(length(t),1);
%noise(15:99,1)=0.1*rand(85,1);
dist=zeros(length(t),1);
%dist(50:69,1)=ones(20,1);

%--------------------------Step----------------------------------------
r =ones(length(t),1);
%-------------------------------------------------------------------------
%...................sine..............................
 %[r,t1]= gensig('sine',length(t)*Ts/2,length(t)*Ts,Ts);
 %-------------------square wave------------------------------------
 %r=[0*ones(floor(length(t)/5),1); 1*ones(floor(length(t)/5),1); 0*ones(floor(length(t)/5),1); -1*ones(floor(length(t)/5),1); 0*ones(length(t)-4*floor(length(t)/5),1)];
  %...................................Step with various jump...................................
  %r=[1*ones(floor(length(t)/6),1); 3*ones(floor(length(t)/6),1); 5*ones(floor(length(t)/6),1); -1*ones(floor(length(t)/6),1); 0*ones(floor(length(t)/6),1); 4*ones(length(t)-5*floor(length(t)/6),1)];

for i=1:length(t)-1
   
     %-----------------------------------------------------------------
        x01pa=x01;
        x02pa=x02;
        U1_(:,i+1)=(U1(1,i)+100)*ones(P,1);
        U2_(:,i+1)=(U2(1,i)+100)*ones(P,1);
        x011pa=x01;
        x021pa=x02;
        x012pa=x01;
        x022pa=x02;
        x011per=x01;
        x021per=x02;
        x012per=x01;
        x022per=x02;
        U1per(:,i+1)=(U1(1,i)+100)*(1+delta)*ones(P,1);
        U2per(:,i+1)=(U2(1,i)+100)*(1+delta)*ones(P,1);
        for j=1:P
            u1pa=U1_(j,i+1);
            u2pa=U2_(j,i+1);
            u11pa=U1_(j,i+1);
            u22pa=U2_(j,i+1);
            u11per=U1per(j,i+1);
            u22per=U2per(j,i+1);
            sim('Modelpa');
            x01pa=x1pa(end);
            x02pa=x2pa(end);
            Y_past(j,i+1)=x2pa(end)-441.2;
            x011pa=x1pa1(end);
            x021pa=x2pa1(end);
            Y_past1(j,i+1)=x2pa1(end);
            x012pa=x1pa2(end);
            x022pa=x2pa2(end);
            Y_past2(j,i+1)=x2pa2(end);
            x012per=x1per2(end);
            x022per=x2per2(end);
            Y_per2(j,i+1)=x2per2(end);
            x011per=x1per1(end);
            x021per=x2per1(end);
            Y_per1(j,i+1)=x2per1(end);
        end
     
     %----------------------------------------------------------------
     
    g1(1:P,i+1)=(Y_per1(1:P,i+1)-Y_past1(1:P,i+1))/(delta);
    b1 = zeros(1,P); b1(1,1)= g1(1,i+1);
    a1 = g1(1:P,i+1);
    G1 = toeplitz(a1,b1);
    G1(:,M) = G1(:,M:P)*ones(P-M+1,1);
    G1 = G1(:,1:M);

    g2(1:P,i+1)=(Y_per2(1:P,i+1)-Y_past2(1:P,i+1))/(delta);
    b2 = zeros(1,P); b2(1,1)= g2(1,i+1);
    a2 = g2(1:P,i+1);
    G2 = toeplitz(a2,b2);
    G2(:,M) = G2(:,M:P)*ones(P-M+1,1);
    G2 = G2(:,1:M);
    G=[G1 G2];
   %------------------------------------------------------------------- 
   
    for j=1:P
        Y_d(j,i+1)=(alpha^j)*y+(1-(alpha)^j)*r(i+1); % Programmed
    end
   
    
    gain_DC = g1(end,i+1);
    gain_DC2 = g2(end,i+1);
    Q = eye(P);
    R1 =((gain_DC2/gain_DC)^2)*gamma*gain_DC^2*eye(M);
    R2=gamma*gain_DC2^2*eye(M);
    R=[R1 zeros(M); zeros(M) R2];

    K=(G'*Q*G+R)\(G'*Q);
    
    x01nm=x01;
    x02nm=x02;
   
    while(1)
        
        Dext(:,i+1)=dext(i+1)*ones(P,1);
        E(:,i+1)=Y_d(:,i+1)-Y_past(:,i+1)-Dext(:,i+1)-Dnl(:,i+1);
        H = 2*(G'*Q*G+R);
        H=(H+H')/2;
        f = -(2*E(:,i+1)'*Q*G)';
        Tl = ones(2*M);
        T = tril(Tl);
        A = [T;-T];
        ub1=up1*ones(M,1)-U1(1,i);
        ub2=up2*ones(M,1)-U2(1,i);
        ub=[ub1; ub2];
        lb1=-low1*ones(M,1)+U1(1,i);
        lb2=-low2*ones(M,1)+U2(1,i);
        lb=[lb1; lb2];
        B=[ub; lb];
        opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');
        dU(:,i+1)= quadprog(H,f,A,B,[],[],[],[],dU01(:,i),opts);
        dU01(:,i+1)=dU(:,i+1);
        dU1(:,i+1)=dU(1:M,i+1);
        dU2(:,i+1)=dU(M+1:2*M,i+1);
        U1(:,i+1)=dU1(:,i+1)+U1(:,i);
        U2(:,i+1)=dU2(:,i+1)+U2(:,i);
        
        U1nm(1:M,i+1)=U1(:,i+1);
        U1nm(M+1:P,i+1)=U1(M,i+1);
        U2nm(1:M,i+1)=U2(:,i+1);
        U2nm(M+1:P,i+1)=U2(M,i+1);
        for j=1:P
            u1nm=U1nm(j,i+1)+100;
            u2nm=U2nm(j,i+1)+100;
            sim('Modelnm');
            x01nm=x1nm(end);
            x02nm=x2nm(end);
            Y_nm(j,i+1)=x2nm(end)-441.2;
        end
         
        Y_m(:,i+1)=G*dU(:,i+1)+Y_past(:,i+1);
        Y_eli(:,i+1)=Y_m(:,i+1)+Dnl(:,i+1);
        old=Dnl(:,i+1);
        
        Dnl(:,i+1)=(1-betta)*Dnl(:,i+1)+betta*(Y_nm(:,i+1)-Y_m(:,i+1));
        n(i+1)=n(i+1)+1;
        if (norm(Dnl(:,i+1)-old))/norm(old)<=0.02 || n(i+1)>=20
            break
        end
        
    end
    
    u1=U1(1,i+1)+100;
    u2=U2(1,i+1)+100;
    sim('Model');
    
    yp=y(end)+noise(i,1)+dist(i,1);
    dext(i+2)=yp-Y_nm(1,i+1)-441.2;
    y=yp-441.2;
    x01=x1(end);
    x02=x2(end);
    yre=[yre; yp];
    yeli=[yeli; Y_eli(1,i+1)];
    ym=[ym; Y_m(1,i+1)];
    u_1=[u_1; u1];
    u_2=[u_2; u2];
    x_1=[x_1; x1(end)];
    
end

figure(1);
subplot(2,2,1:2);
plot(yre,'b');
hold on
plot(r+441.2,'r');
grid on
legend('y','r');
title('Response of the system');
xlabel('sample');
%------------------------
subplot(2,2,2);
plot(yre,'b');
hold on
plot(yeli+441.16,'r');
hold on
plot(ym+441.2,'m');
grid on   
legend('y','r');
title('Yeli & Ym with bias and Yp');
legend('Yprocess','Yextended linear model','Ymodel');
xlabel('sample');
%------------------------
subplot(2,2,3);
plot(u_1,'b');
grid on
xlabel('sample');
title('Control law for input 1');
%-------------------------
subplot(2,2,4);
plot(u_2,'b');
grid on
xlabel('sample');
title('Control law for input 2');
%--------------------------
% subplot(3,2,5:6);
% plot(x_1)
% grid on
% xlabel('sample');
% title('state one');
% subplot(3,2,4);
% plot(dist,'b');
% grid on
% xlabel('sample');
% title('disturbance');
% plot(noise,'b');
% grid on
% xlabel('sample');
% title('noise');
% figure
% plot(x_1);

    
    
    
    
        
    
