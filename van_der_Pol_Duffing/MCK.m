function [M,C,K,N_dof]=MCK
%% 说明
  %% 该函数功能的简单介绍：
  %输出M，C，K，N_dof。便于我们修改模型
  %% 该函数实现的方法介绍
  %把论文里头的各种参数输入进来就好
%% 实现
alpha=0.2;w=0.5;
M=1;
C=-alpha;
K=w^2;
N_dof=size(M,1);
