%% 说明
  %% 该程序功能的简单介绍：
  %目的是验证Jacibo是正确的
  %% 该程序实现的方法介绍
  %很简单不完全，还没写完
%% 实现
clc
clear
M=[1 0 0;0 1 0;0 0 1];C=[1 0 0; 0 0 0; 0 0 1];K=[1 0 0; 0 2 0;0 0 3];
N_dof=size(M,1);
n=2;
for j=1:2*N_dof
    for i=1:n
        da=zeros(n,2*N_dof);
        da(i,j)=1; da_1=matrixder(da,1);da_2=matrixder(da,2);
        part_1=matri_mul_ma(M,da_2); solution=part_1;
        solution=solution(1:n,:);
        solution=matrixtra(solution);
        solution=arrange_column(solution);
       S(:,i+(j-1)*n)=solution;
    end
end
y=zeros((2*n-1)*N_dof,(2*n-1)*N_dof);
for i=1:N_dof
    y(:,(i-1)*(2*n-1)+1:(i-1)*(2*n-1)+n)=S(:,(i-1)*(2*n)+1:(i-1)*(2*n)+n);
    y(:,(i-1)*(2*n-1)+n+1:(i-1)*(2*n-1)+2*n-1)=S(:,(i-1)*(2*n)+n+2:(i-1)*(2*n)+2*n);
end
