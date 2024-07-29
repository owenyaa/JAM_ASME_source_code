function y=matrixsum_2(A,B)
%% 说明
  %% 该函数功能的简单介绍：
  %是求三角函数系数矩阵（cosnt,sinnt前面的系数）和的函数
  %A和B是需要求和系数矩阵，也是n*2的矩阵
  %% 该函数实现的方法介绍
  %先判断A、B两矩阵的长度
  %在新建一个矩阵，长度取A,B长度的最大值
  %取长度大的矩阵为A，小的为B
  %新的矩阵里头的每一个元素都等于AB的元素之和

%% 实现
  %求解矩阵的维度
  L1=size(A,1);L2=size(B,1);
  L3=size(A,2)/2;
  L=max(L1,L2);
  solution=zeros(L,2*L3);  
  
  %交换矩阵
  if L1 < L2
      A1=B;
      B=A;
      A=A1;
      L2=L1;
      L1=L;
  end 
  
  %%求解矩阵的元素
  for i=1:L2
      for j=1:2*L3
          solution(i,j)=A(i,j)+B(i,j);
      end
  end
  for i=L2+1:L1
      for j=1:2*L3
          solution(i,j)=A(i,j);
      end
  end
 for i=1:L3
 solution(1,2*i)=0;
 end
  y=solution;
  %% 以下是可能会用得到的代码
  %A=[1 2 3 4 5 6 7 8;2 2 3 4 5 6 7 8 ;3 2 3 4 5 6 7 8;4 2 3 4 5 6 7 8]
  %B=[1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1]
  %y=matrixsum_2(A,B)