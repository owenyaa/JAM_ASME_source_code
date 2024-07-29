%% 说明
  %% 该程序功能的简单介绍：
  %用IHB法求解非线性方程组（含外激励的周期振动问题）
  %% 该程序实现的方法介绍
  %给定初值，通过循环求解增量的方法迭代求解方程的自激运动及其频率
  %具体理论方法参见IHB法，这里面只给出计算机实现过程中的重点
%% 实现
clear;%clc;
close all;
tic
global N_harm Q N_dof
%w0是自激振动的频率也是未知数，a是三角级数的系数也是未知数，n是所用的谐波阶数,Q是风速（不一定用到），N_dof是自由度个数
N_harm=20;N_dof=3;
% QQ=2.3:0.05:4.5;
% QQ=2.3:0.01:3.4;%P3
% QQ=4.5:-0.01:2.4;%P1
QQ=3:-0.01:2.4;%P2_1
% QQ=3:0.05:3.2;%P2_2

parameter_a=zeros(N_harm,2*N_dof);
% % a(2,3)=-2;a(2,4)=1.6;%P3
% % a(2,3)=0;a(2,4)=-4.6;%P1
% a(2,3)=-2.4;a(2,4)=-3.4;%P2
% % a(2,:)=[-0.751428351337642,-0.829293669552665,-0.124060495020285,-0.164596039902716,-0.167687644847091,-0.200533944785485];%P2
% parameter_a(2,:)=[-0.137716994609563   0.081667657228171  -0.051267163033702   0.026281681668196  -0.062209577472483   0.033821920088153];%P3 Q=2.3
% parameter_a(2,:)=[-2.529278464743486   1.811833020583184  -0.455319079257941   0.285457945854791  -0.779600182786925   0.528170101402331];%P1 Q=4.5
parameter_a(2,:)=[-0.971936911179476  -0.756650814395171  -0.179366762632965  -0.163605789072816  -0.249976193884288  -0.209328091146364];%P2 Q=3
% w0=0.4;
% w0=0.450731040778181;%P3 Q=2.3
% w0=0.450609296757463;%P1 Q=4.5
w0=0.425949568961625;%P2 Q=3
for ii=1:length(QQ)
    Q=QQ(ii)
    [M,C,K,N_dof]=MCK(Q);
    er=1;tol=1e-7;iter=0;
    while(er>=tol&&iter<1000)
        part=IHB_part(parameter_a,w0);
        R=part(1).vector; Ja=part(2).vector; Jw=part(3).vector;
        er=norm(R);
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
        parameter_a=parameter_a+da;
        iter=iter+1;
    end
    A_11=sqrt(parameter_a(2,1)^2+parameter_a(2,2)^2);
    A_21=sqrt(parameter_a(2,3)^2+parameter_a(2,4)^2);
    A_31=sqrt(parameter_a(2,5)^2+parameter_a(2,6)^2);
    every(ii).w=w0;every(ii).Q=Q;
    every(ii).parameter_a=parameter_a;
    every(ii).A_11=A_11;every(ii).A_21=A_21;
    every(ii).A_31=A_31;
end
% A_21=sqrt(a(2,3)^2+a(2,4)^2)