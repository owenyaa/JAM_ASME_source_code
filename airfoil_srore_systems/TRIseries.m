function y=TRIseries(n)
%% 说明
  %% 该函数功能的简单介绍：
  %用来求解n阶的三角矩阵的函数
  %输出的是列向量
  %% 该函数实现的方法介绍
  %【a0,0;cos(t),sin(t);cos(2*t);sin(2*t).....】
%% 实现
syms t
S=sym(zeros(n,2));
for i=1:n
    if i==1
        S(i,1)=1;
    else
        S(i,1)=cos((i-1)*t);
    end
end
for i=1:n
    if i==1
        S(i,2)=0;
    else
        S(i,2)=sin((i-1)*t);
    end
end
y=sym(zeros(2*n-1,1));
y(1:n,1)=S(:,1);
y(n+1:end,1)=S(2:end,2);

