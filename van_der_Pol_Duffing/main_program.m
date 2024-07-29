%% 说明
%% 该程序功能的简单介绍：
%用IHB法求解非线性方程组（含外激励的周期振动问题）
%% 该程序实现的方法介绍：
%给定初值，通过循环求解增量的方法迭代求解方程的运动及其频率
%具体理论方法参见IHB法，这里面只给出计算机实现过程中的重点
%% 实现
clear;clc;close all;
tic
global n M C K N_dof
%w0是自激振动的频率也是未知数，a是三角级数的系数也是未知数，n是所用的谐波阶数,Q是风速（不一定用到），N_dof是自由度个数
n=20;record=[];istep=1;

%% 系统参数
alpha=0.2;w=0.5;
M=1;
C=-alpha;
K=w^2;
N_dof=size(M,1);

w0=0.45;
% load a_w_0.45.mat
% [a,~]=get_initial(n,w0);
% load a_initial.mat
a=zeros(n,N_dof*2);
% a(2,:)=[a1,a2];
a(2,:)=[1.2,0.2];
er=1;tol=1e-10;
dw=0;iter=0;
while(er>=tol&&iter<1001)
    part=IHB_part(a,w0);
    R=part(1).vector; Ja=part(2).vector;
    %% TIHB method
    [u,s,v]=csvd(Ja);
    lambda=l_curve(u,s,R);
    da=tikhonov(u,s,v,R,lambda);
    
    %% the traditional IHB medhod
    %     da=Ja\R;
    da=arrange_column_inv(da);
    da=matrixtra_inv(da);
    er=norm(R)
    w0=w0+dw;
    a=a+da;
    iter=iter+1;
    if er>1&&iter>500
        '不收敛'
        break
    end
end
% a
sqrt(a(2,1)^2+a(2,2)^2)
toc

dt=0.01;
Tdata=0:dt:40;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
Harm_parameter_a=a(2:end,:);
for j=1:N_dof
    for i=1:n-1
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(i)*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*(i)*Harm_parameter_a(i,2*j)*cos((i)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(i))^2*Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)-(w0*(i))^2*Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
    end
    x(j,:)=x(j,:)+a(2*j-1,1);
end
beta=0.3;f=1;r=0.2;
F=beta*x.^3+r*x.^2.*dx;
f_ex=f*cos(w0*Tdata);
% residual1(1:N_dof,:)=M*ddx+C*dx+K*x;
residual(1:N_dof,:)=M*ddx+C*dx+K*x+F-f_ex;
figure;
plot(Tdata,residual(1,:),'b--','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
%
% % dt=tt(2)-tt(1);
% data=residual(1,:);
% N=length(data);
% %tt=0:dt:(N-1)*dt;tt=tt';
% N_fft=2^16;
% Y=fft(data,N_fft);
% Pyy=Y.*conj(Y)/N_fft;
% f=2*pi/dt*(0:N_fft/2)/N_fft;
% figure
% plot(f,(Pyy(1:(N_fft/2+1))))
% title('Frequency content of y')
% xlabel('frequency (Hz)')


figure;
plot(Tdata,x,'r-','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% %
%% 数值解法
Tend=1000;
dt=0.01;
Tspan=0:dt:Tend;
% y0=ones(2*N_dof,1);
% y0=[1,1];
y0=[x(1,1),dx(1,1)];
options=odeset('RelTol',1e-10);
[T,X]=ode45(@(t,x) odefun(w0,t,x), Tspan, y0, options);
figure;
plot(Tdata,x,'r-','LineWidth',1.5);
hold on;
plot(Tspan(1:40:end),X(1:40:end,1),'k.','MarkerSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

figure;
plot(x(1,1:period*100),dx(1,1:period*100),'r-','LineWidth',1.5);
hold on;
plot(X(1:20:period*100,1),X(1:20:period*100,1+N_dof),'k.','MarkerSize',15)
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

%% ODE函数
function y=odefun(w0,t,x)
[M,C,K,N_dof]=MCK;
beta=0.3;f=1;r=0.2;
F=beta*x(1)^3+r*x(1)^2*x(2);
f_ex=f*cos(w0*t);
y=zeros(2*N_dof, 1);
y(1:N_dof)=x(N_dof+1:end);
y(N_dof+1:end)=M\(f_ex-K*x(1:N_dof)-C*x(N_dof+1:end)-F);
end
