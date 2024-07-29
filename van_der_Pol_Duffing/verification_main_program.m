%% 说明
  %% 该程序功能的简单介绍：
  %用IHB法求解非线性方程组（含外激励的周期振动问题）
  %% 该程序实现的方法介绍
  %给定初值，通过循环求解增量的方法迭代求解方程的运动及其频率
  %具体理论方法参见IHB法，这里面只给出计算机实现过程中的重点
%% 实现
clear
clc
close all
tic
global w0 a n M C K N_dof
%w0是自激振动的频率也是未知数，a是三角级数的系数也是未知数，n是所用的谐波阶数,Q是风速（不一定用到），N_dof是自由度个数
n=20;iter=0;
[M,C,K,N_dof]=MCK;
w0=0.3;
% [a,~]=get_initial(n);a(5:end,:)=zeros(n-5+1,2*N_dof);

a=zeros(n,2*N_dof);
a(1,:)=[0,0];
a1=0.3;a2=-0.9;
a(2,:)=[a1,a2];

%% 循环求解增量不断接近精确解
er=1;tol=1e-6;
dw=0;
while(er>=tol&&iter<1001)
    part=IHB_part(M,C,K);
    R=part(1).vector; Ja=part(2).vector; Jw=part(3).vector;
    er=norm(R)
    da=Ja\R;
    %     [u,s,v]=csvd(Ja);
    %     lambda=l_curve(u,s,R);
    %     da=tikhonov(u,s,v,R,lambda);
    da=arrange_column_inv(da);
    da=matrixtra_inv(da);
    w0=w0+dw;
    a=a+da;
    iter=iter+1;
end
toc