function y=matrixtra(A)
%% 说明
  %% 该函数功能的简单介绍：
  %将三角矩阵变换为多列向量的形式
  %输入的是A（n*2N_dof),n是谐波数，N_dof是自由度数
  %% 该函数实现的方法介绍
  %第一列不变，第二列从第二行开始放在第一列下方
  %有多少个自由度就有多少列
%% 实现
[L1,L2]=size(A);L2=L2/2;
solution=zeros(2*L1-1,L2);
for i=1:L2
    solution(1:L1,i)=A(:,2*i-1);
    solution(L1+1:end,i)=A(2:end,2*i);
end
y=solution;
%% 以下是可能用的到的代码
%A=[1 2 3 4 5 6 7 8;2 2 3 4 5 6 7 8 ;3 2 3 4 5 6 7 8;4 2 3 4 5 6 7 8]