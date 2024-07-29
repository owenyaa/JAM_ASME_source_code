function part=IHB_part(a,w0)
%% 说明
  %% 该函数功能的简单介绍：
  %只能用来求解自激振动
  %计算出模型残差，Ja，Jw的三角矩阵
  %M,C,K是质量阻尼刚度矩阵，F是非线性项(要自己提前算出来带进来）
  %输出残差是列向量，Ja是矩阵，Jw是列向量，是arrange_column的输出值
  %% 该函数实现的方法介绍
  %将所有的参数带入方程中求解，然后输出
%% 实现
global  n M C K N_dof

da=zeros(n,2*N_dof);
a_1=matrixder(a,1);a_2=matrixder(a,2);
F=nonterm(w0,a,da);
%% 残差
R(1).matrix=matri_mul_ma(M,-w0^2*a_2); R(2).matrix=matri_mul_ma(C,-w0*a_1); R(3).matrix=matri_mul_ma(K,-a); R(4).matrix=-F(1).part;
R_matrix=matrixsum(R);
R_vector=matrixtra(R_matrix(1:n,:));
R_vector=arrange_column(R_vector);
%% Ja
Ja=jacobi(a,w0);
%% Jw
Jw(1).matrix=matri_mul_ma(M,2*w0*a_2); Jw(2).matrix=matri_mul_ma(C,a_1); Jw(3).matrix=F(3).part;
Jw=matrixsum(Jw); 
Jw=matrixtra(Jw(1:n,:));
Jw=arrange_column(Jw);
%% 组合
part(1).vector=R_vector; part(2).vector=Ja; part(3).vector=Jw;