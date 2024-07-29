function y=matrixder(A,n)
%% 说明
  %% 该函数功能的简单介绍：
  %用来对三角函数系数矩阵的求导（cosnt,sinnt前面的系数）
  %A是需要求导的系数矩阵,n是求导的次数
  %输入和输出的都是（谐波数*2倍自由度）维度的矩阵
  %% 该函数实现的方法介绍
  %三角函数求导会有额外的系数，求一次导就会多出一个n，且cos变-sin，sin变cos
  %通过n/2的余数来判断求了奇次导还是偶次导
  %偶数列的性质相同（与第二列相同），奇数列的性质相同（与第一列相同）
 %% 实现
 %循环求解解矩阵系数
 [L1,L2]=size(A);
 solution=zeros(L1,L2);
 remainder=mod(n,2); 
 for i=1:L1
     for j=1:L2
         if i==1 && n>=1
             solution(i,j)=0;
         end
         if remainder==1
             solution(i,j)=(i-1)^n*(-1)^((n-3+2*(2-mod(j,2)))/2)*A(i,j+(-1)^(mod(j,2)+1));
         else
             solution(i,j)=(i-1)^n*(-1)^(n/2)*A(i,j);
         end
     end
 end
 for i=1:L2/2
 solution(1,2*i)=0;
 end
 y=solution;