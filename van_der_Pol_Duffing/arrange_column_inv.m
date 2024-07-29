function y=arrange_column_inv(A)
%% 说明
  %% 该函数功能的简单介绍：
  %该函数将一列向量转换为多列矩阵（并不是三角矩阵形式）
  %输入列向量，输出多列矩阵
  %% 该函数实现的方法介绍
  %根据自由度平均列向量
  global N_dof
  L1=size(A,1);
  n=L1/N_dof;
  solution=zeros(n,N_dof);
  for i=1:N_dof
      solution(1:end,i)=A(1+(i-1)*n:i*n,:);
  end
  y=solution;