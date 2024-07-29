function F=nonterm(w0,a,da)
%% 说明
  %% 该函数功能的简单介绍：
  %这是要求解的自激振动方程中的非线性项，如果要换模型，需要改变这里面的值
  %w0和a是全局变量，da是增量(a,和da都是三角矩阵的形式）
  %输出的数组每一个元素都是多列向量，可能不同维所以用数组表示（同matrixtra函数输出的值一致）
  %% 该函数实现的方法介绍
  %F被分成了三个部分（1、残差 2、Ja 3、Jw）
%% 实现
global N_dof n
beta=0.3;r=0.2;f=1;

%将a拆成三个单自由度去运算
for i=1:N_dof
    X(i).part=a(:,2*i-1:2*i);
    dX(i).part=da(:,2*i-1:2*i);
end
for i=1:N_dof
    X_1(i).part=matrixder(X(i).part,1);
    X_2(i).part=matrixder(X(i).part,2);
    dX_1(i).part=matrixder(dX(i).part,1);
    dX_2(i).part=matrixder(dX(i).part,2);
end

%% 残差
term1_1(1).matrix=beta*X(1).part; term1_1(2).matrix=X(1).part; term1_1(3).matrix=X(1).part; Term1(1).matrix=matrixmul(term1_1);
term1_2(1).matrix=w0*r*X(1).part; term1_2(2).matrix=X(1).part; term1_2(3).matrix=X_1(1).part; Term1(2).matrix=matrixmul(term1_2);
term1_3=zeros(n,2); term1_3(2,1)=f; Term1(3).matrix=-term1_3;%外激励
Term1=matrixsum(Term1);
F(1).part=Term1;
% F(3).part=zeros(n,2*N_dof);

%% Ja
term2_1(1).matrix=3*beta*X(1).part; term2_1(2).matrix=X(1).part; term2_1(3).matrix=dX(1).part; term2_1=matrixmul(term2_1); Term2(1).matrix=term2_1;
term2_2(1).matrix=w0*r*X(1).part; term2_2(2).matrix=X(1).part; term2_2(3).matrix=dX_1(1).part; term2_2=matrixmul(term2_2); Term2(2).matrix=term2_2;
term2_3(1).matrix=2*w0*r*X(1).part; term2_3(2).matrix=X_1(1).part; term2_3(3).matrix=dX(1).part; term2_3=matrixmul(term2_3); Term2(3).matrix=term2_3;
F(2).part=matrixsum(Term2);

%% Jw
term3(1).matrix=r*X(1).part; term3(2).matrix=X_1(1).part; term3(3).matrix=X_1(1).part; 
Term3=matrixmul(term3); 
F(3).part=Term3;
