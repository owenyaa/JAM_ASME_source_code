% clc
% clear
% close all
global Q N_harm N_dof
N_harm=20;N_dof=3;
% Q0=2.4;dQ=0.1;
% record=zeros(5,1);


Q=3;
a1=10*(2*rand(1)-1);a2=10*(2*rand(1)-1); a3=10*(2*rand(1)-1);a4=10*(2*rand(1)-1); a5=10*(2*rand(1)-1);a6=10*(2*rand(1)-1);
% a1=-7.0892;a2=-7.2786;a3=7.3858;a4=1.5941;a5=0.99972;a6=-7.1009;
% a1=-8.0709;a2=-7.3605;a3=8.8410;a4=9.1227;a5=1.5042;a6=-8.8044;
% a1=0.6847;a2=0.1604;a3=0.3036;a4=0.3222;a5=0.1345;a6=0.0734;%这是一个很奇怪的解
y0=[a1,a2,a3,a4,a5,a6];

Tend=50000;
dt=0.01;
Tspan=0:dt:Tend;


options=odeset('RelTol',1e-10);
[T,X]=ode45(@(t,x) odefun(x), Tspan, y0, options);

% load Q_2.4_11.3_numerical.mat
% Q=11.3;dQ=0.1;
% X=record_x(end).part;
% y0=X(end,:);
% for i=91:97
%     Q=Q+dQ;
% % y0=ones(2*N_dof,1);
% Tend=5000;
% dt=0.01;
% Tspan=0:dt:Tend;
%
%
% options=odeset('RelTol',1e-10);
% [T,X]=ode45(@(t,x) odefun(x), Tspan, y0, options);
% record(i,1)=Q;
% record_x(i).part=X;
% i
% y0=X(end,:);
% end

% figure;
% plot(X(end-50/dt:end,1),X(end-50/dt:end,1+N_dof),'r.','LineWidth',1)
% hold on
plot(X(end-500/dt:end,2),X(end-500/dt:end,2+N_dof),'r.','LineWidth',1)
% hold on
% plot(X(end-50/dt:end,3),X(end-50/dt:end,3+N_dof),'r.','LineWidth',1)




function y=odefun(x)
global Q
k_h3=0;k_alpha3=10;k_beta3=0;
[M,C,K,N_dof]=MCK(Q);
F=[k_h3*x(1)^3;k_alpha3*x(2)^3;k_beta3*x(3)^3];
y=zeros(2*N_dof, 1);
y(1:N_dof)=x(N_dof+1:end);
y(N_dof+1:end)=M\(-K*x(1:N_dof)-C*x(N_dof+1:end)-F);

end
