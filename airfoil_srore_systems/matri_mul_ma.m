function y=matri_mul_ma(A,B)
%% 说明
  %% 该函数功能的简单介绍：
  %是三角函数系数矩阵和方阵相乘的函数A
  %A是方阵（N_dof*N_dof),B是三角系数矩阵（n*2N_dof)
  %输出的是三角系数矩阵（n*2N_dof)
  %% 该函数实现的方法介绍
  %
%% 实现
[n,N_dof]=size(B);
N_dof=N_dof/2;
for i=1:N_dof
    X(i).part=B(:,2*i-1:2*i);
end
for i=1:N_dof
    solution(i).part=A(i,1)*X(1).part;
end
for i=1:N_dof
    for j=1:N_dof-1
        solution(i).part=matrixsum_2(solution(i).part,A(i,j+1)*X(j+1).part);
    end
end
y=zeros(n,2*N_dof);
for i=1:N_dof
    y(:,2*i-1:2*i)=solution(i).part;
end