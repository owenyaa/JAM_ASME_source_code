function y=matrixtra_inv(A)
%% 说明
  %% 该函数功能的简单介绍：
  %将多列向量变换为三角系数矩阵
  %% 该函数实现的方法介绍
  %取列向量长度的一半（要向上取整），第二列要添一个0
%% 实现

[L1,L2]=size(A);
L1=ceil(L1/2);
solution=zeros(L1,2*L2);
for i=1:L2
    solution(:,2*i-1)=A(1:L1,i);
    solution(2:end,2*i)=A(L1+1:end,i);
end
y=solution;
%% 以下是可能用的到的代码
% y=matrixtra_inv(matrixtra(A))