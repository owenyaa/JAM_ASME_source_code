%% ˵��
  %% �ó����ܵļ򵥽��ܣ�
  %��IHB���������Է����飨���⼤�������������⣩
  %% �ó���ʵ�ֵķ�������
  %������ֵ��ͨ��ѭ����������ķ���������ⷽ�̵��Լ��˶�����Ƶ��
  %�������۷����μ�IHB����������ֻ���������ʵ�ֹ����е��ص�
%% ʵ��
clear;clc;%close all;
tic
global N_harm Q N_dof
%w0���Լ��񶯵�Ƶ��Ҳ��δ֪����a�����Ǽ�����ϵ��Ҳ��δ֪����n�����õ�г������,Q�Ƿ��٣���һ���õ�����N_dof�����ɶȸ���
N_harm=20;Q=3;

[M,C,K,N_dof]=MCK(Q);
parameter_a=zeros(N_harm,2*N_dof);

% parameter_a(2,3)=10*(2*rand(1)-1);parameter_a(2,4)=10*(2*rand(1)-1);
%TIHB
% parameter_a(2,3)=1.2;parameter_a(2,4)=3;%P1 ��
% parameter_a(2,3)=1;parameter_a(2,4)=1.2;%P2 ��
% parameter_a(2,3)=-1.2;parameter_a(2,4)=4.2;%P3 ��

%IHB
% parameter_a(2,3)=1.2;parameter_a(2,4)=0.2;%P1 ��
% parameter_a(2,3)=3.2;parameter_a(2,4)=7.8;%P2 ��
parameter_a(2,3)=4.4;parameter_a(2,4)=4.4;%P3 ��

% load a_Q_3_2.mat;
w0=0.4;
% w0=0.559031674222567;
% load Q_3.2_monocyclic.mat
% [parameter_a,w0]=get_initial(N_harm);
% a(4:end,:)=0;
% %% ѭ������������Ͻӽ���ȷ��
er=1;tol=1e-7;iter=0;
tic
while(er>=tol&&iter<1000)
    part=IHB_part(parameter_a,w0);
    R=part(1).vector; Ja=part(2).vector; Jw=part(3).vector;
    er=norm(R)
    %     Ja(:,2)=Jw;
    Ja(:,4)=Jw;
    [u,s,v]=csvd(Ja);
    lambda=l_curve(u,s,R);
    da=tikhonov(u,s,v,R,lambda);
    da=Ja\R;
    dw=da(4);
    da(4)=0;
    
    
    %     Ja(:,2)=Jw;
    %     %Ja(:,4)=Jw;
    %     %     [u,s,v]=csvd(Ja);
    %     %     lambda=l_curve(u,s,R);
    %     %     da=tikhonov(u,s,v,R,lambda);
    %     da=Ja\R;
    %     dw=da(2);
    %     da(2)=0;
    
    
    da=arrange_column_inv(da);
    da=matrixtra_inv(da);
    w0=w0+dw;
    parameter_a=parameter_a+da;
    iter=iter+1;
end
aa=toc;
% norm(a)
w0=abs(w0);
A_11=sqrt(parameter_a(2,3)^2+parameter_a(2,4)^2)
% A_21=sqrt(a(2,3)^2+a(2,4)^2)
toc
%% ��ͼ
dt=0.01;
Tdata=0:dt:1000;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm-1
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(i)*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*(i)*Harm_parameter_a(i,2*j)*cos((i)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(i))^2*Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)-(w0*(i))^2*Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
    end
    x(j,:)=x(j,:)+parameter_a(2*j-1,1);
end

k_h3=0;k_alpha3=10;k_beta3=0;
F=[k_h3*x(1,:).^3;k_alpha3*x(2,:).^3;k_beta3*x(3,:).^3];
residual=M*ddx+C*dx+K*x+F;

figure;
plot(Tdata,x(2,:),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);



% axis([1990 2000 -inf inf])

figure;
plot(x(2,end-floor(2*pi/w0/dt):1:end),dx(2,end-floor(2*pi/w0/dt):1:end),'k-','LineWidth',1.5);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
hold on
% plot(x(3,:),dx(3,:),'k','LineWidth',1.5);
% hold on
h1=legend('$$a_1^2$$','$$b_1^2$$');%,);
% h1=legend('$$c_1$$','$$s_1$$','$$\beta$$');
set(h1,'Interpreter','latex','FontSize',15);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% 

%% ��ֵ�ⷨ
% Tend=20000;
% dt=0.01;
% Tspan=0:dt:Tend;
% period=2*pi/w0;
% n_dot=floor(period/dt);
% % y0=ones(2*N_dof,1);
% y0=[x(1,end);x(2,end);x(3,end);dx(1,end);dx(2,end);dx(3,end)];
% options=odeset('RelTol',1e-8,'AbsTol',1e-8);
% [T,X]=ode45(@(t,x) odefun(x), Tspan, y0, options);
% figure;
% % plot(X(end-n_dot:20:end,1),X(end-n_dot:20:end,1+N_dof),'k.','MarkerSize',12)
% % hold on
% plot(X(end-n_dot:35:end,2),X(end-n_dot:35:end,2+N_dof),'k.','MarkerSize',12)
% set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% % figure;
% % plot(Tspan,X(:,1),'r-','LineWidth',1.5);
% % hold on
% % plot(Tdata,x(1,:),'k-','LineWidth',1.5);
% % set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
% 
% % hold on
% % plot(X(end-n_dot:30:end,3),X(end-n_dot:30:end,3+N_dof),'k.','MarkerSize',12)
% %% ODE����
% function y=odefun(x)
% global Q
% % alpha_f=0/180*pi;deta=1/180*pi;M_0=0.5/180*pi;
% k_h3=0;k_alpha3=10;k_beta3=0;
% [M,C,K,N_dof]=MCK(Q);
% F=[k_h3*x(1)^3;k_alpha3*x(2)^3;k_beta3*x(3)^3];
% y=zeros(2*N_dof, 1);
% y(1:N_dof)=x(N_dof+1:end);
% y(N_dof+1:end)=M\(-K*x(1:N_dof)-C*x(N_dof+1:end)-F);
% 
% end