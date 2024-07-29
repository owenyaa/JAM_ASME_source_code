%% 说明
  %% 该程序功能的简单介绍：
  %用IHB法求解非线性方程组（自激问题）
  %% 该程序实现的方法介绍
  %给定初值，通过循环求解增量的方法迭代求解方程的自激运动及其频率
  %具体理论方法参见IHB法，这里面只给出计算机实现过程中的重点
%% 实现
clear
clc
close all
tic
global w0 a n M C K N_dof Q t_a
%w0是自激振动的频率也是未知数，a是三角级数的系数也是未知数，n是所用的谐波阶数,Q是风速（不一定用到），N_dof是自由度个数
n=20;Q=6.2851*0.75;
[M,C,K,N_dof]=MCK;
% w0=0.079611595418043;
[a,w0]=get_initial(n);
% n_1=floor(n/2);
% for i=1:n_1
%     a(2*i-1,:)=zeros(1,2*N_dof);
% end
% a=zeros(n,2*N_dof);
% a(2,:)=[10   10  10  10 10  -10];
% w0=1;
% a(4,:)=[0.000096381114544    0.000099215253183  -0.002752690285901  0.002906575466455  -0.023458216001988  -0.036994651364661];
% a(6,:)=[-0.000069699852260   0.000078790021580   0.000318261918995 -0.000498961353908   0.001513551073024   0.002691599690736];
% a(8,:)=[0.000013460056720   -0.000035264751157  -0.000043621573600  0.000165728562619  -0.000349605006046  -0.000388071403407];
% a(10,:)=[-0.000000132412546   0.000014790167236  -0.000004786841998 -0.000062280014010   0.000109131050396   0.000063513633537];

da=[];

er=1;tol=1e-6;iter=0;
%% 循环求解增量不断接近精确解

    
while(er>=tol&&iter<1000)
    part=IHB_part(M,C,K);
    R=part(1).vector; Ja=part(2).vector; Jw=part(3).vector;
    er=norm(R)
    Ja(:,4)=Jw;
    [u,s,v]=csvd(Ja);
    lambda=l_curve(u,s,R);
    da=tikhonov(u,s,v,R,lambda);
%     da=Ja\R;
    dw=da(4);
    da(4)=0;
    da=arrange_column_inv(da);
    da=matrixtra_inv(da);
    w0=w0+dw;
    a=a+da;
    iter=iter+1;
end

iter
toc




Tdata=0:0.01:400;
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
figure;
plot(x(1,:),dx(1,:),'r','LineWidth',1);
% hold on
figure;
plot(x(2,:),dx(2,:),'b','LineWidth',1);
% hold on
% plot(x(3,:),dx(3,:),'k','LineWidth',1);

% figure;
% plot(Tdata,x(1,:),'r-','LineWidth',1);
% hold on;
% plot(Tdata,x(2,:),'b-','LineWidth',1);
% hold on
% plot(Tdata,x(2,:),'k-','LineWidth',1);

