 function y=jacobi(a,w0)
 %% 说明
  %% 该函数功能的简单介绍：
  %求一个自激振动的Jacobi矩阵
  %输出的Jacobi矩阵应该是(2n-1)*（N_dof(2n-1))维
  %% 该函数实现的方法介绍
  %通过求偏导（由于是线性的，可以直接把某一个值令为一，其他的都为零），来求Jacobi矩阵
%% 实现
global n M C K N_dof

for j=1:2*N_dof
    
    for i=1:n
        da=zeros(n,2*N_dof);
        da(i,j)=1; da_1=matrixder(da,1);da_2=matrixder(da,2);
        part_1=matri_mul_ma(M,w0^2*da_2); Ja(1).matrix=part_1;
        part_2=matri_mul_ma(C,w0*da_1); Ja(2).matrix=part_2;
        part_3=matri_mul_ma(K,da); Ja(3).matrix=part_3;
        F=nonterm(w0,a,da); 
        part_4=F(2).part; Ja(4).matrix=part_4;
        solution=matrixsum(Ja);
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
    
