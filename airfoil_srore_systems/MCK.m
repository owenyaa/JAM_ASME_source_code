function [M,C,K,N_dof]=MCK(Q)
%% 说明
  %% 该函数功能的简单介绍：
  %输出M，C，K。便于我们修改模型
  %% 该函数实现的方法介绍
  %把论文里头的各种参数输入进来就好
%% 实现


mu=12.8;
mu_beta=4.0;
x_alpha=0.15;
x_beta=0.2;
L_1=0.18;
r_alpha=sqrt(0.3);
r_beta=sqrt(0.89);
c_h=0.2;
c_alpha=0.2;
c_beta=0;
w_h=34.6;
w_alpha=88;
w_beta=60;
a_con=-0.41;

M=[mu+mu_beta, mu*x_alpha+mu_beta*x_beta-mu_beta*L_1, mu_beta*x_beta;...                                                                                                                                          
    mu*x_alpha+mu_beta*x_beta-mu_beta*L_1, mu*(r_alpha)^2+mu_beta*(r_beta)^2+mu_beta*(L_1)^2-2*mu_beta*x_beta*L_1, mu_beta*(r_beta)^2-mu_beta*x_beta*L_1;...
    mu_beta*x_beta, mu_beta*(r_beta)^2-mu_beta*x_beta*L_1, mu_beta*(r_beta)^2];

C=[c_h, 0, 0; 0, c_alpha, 0; 0, 0, c_beta];

K=[mu*(w_h/w_alpha)^2, 2*Q, 0; 0, mu*r_alpha^2-(1+2*a_con)*Q, 0; 0, 0, mu_beta*(r_beta)^2*(w_beta/w_alpha)^2];
N_dof=size(M,1);
%% 以下是可能会用到的代码
% M=[1,2,3;...
%     2,3,4;...
%     3,4,5];