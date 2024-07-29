function y=matrixsum(every)
%% 说明
  %% 该函数功能的简单介绍：
  %是求三角函数系数矩阵（cosnt,sinnt前面的系数）和的函数
  %every是需要求和系数矩阵数组，里头有多个矩阵，都是n*2的矩阵
  %% 该函数实现的方法介绍
  %先判断every中矩阵的长度
  %在新建一个矩阵，长度取长度的最大值
  %再用两个系数矩阵相加的小函数一直求和

%% 实现
  %求解矩阵的维度
  N=size(every,2);
  len=zeros(N,1);
  for i=1:N
      len(i)=size(every(i).matrix,1);
      L1=max(len);
  end
  L2=size(every(1).matrix,2);
  S=zeros(L1,L2);
  %不断用子函数求解
  for i=1:N
      S=matrixsum_2(S,every(i).matrix);
  end
  
  y=S;