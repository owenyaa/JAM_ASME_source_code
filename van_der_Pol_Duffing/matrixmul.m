function y=matrixmul(every)
%% 说明
  %% 该函数功能的简单介绍：
  %是求多个三角函数系数矩阵（cosnt,sinnt前面的系数）相乘的函数
  %% 该函数实现的方法介绍
  %every中的元素一个一个的相乘
%% 实现
  N=size(every,2);
  S=every(1).matrix;
  for i=2:N
      S=matrixmul_2(S,every(i).matrix);
  end
  y=S;