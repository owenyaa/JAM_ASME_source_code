function y=matrixmul_2(A,B)
%% 说明
  %% 该函数功能的简单介绍：
  %是求三角函数系数矩阵（cosnt,sinnt前面的系数）相乘的函数
  %A和B是需要求乘系数矩阵，也是n*2的矩阵
  %% 该函数实现的方法介绍
  %先判断A、B两矩阵的长度，来确定解矩阵的长度
  %分成四个部分来求解解矩阵中的元素
  
%% 实现
L1=size(A,1);L2=size(B,1);
L3=size(A,2)/2;
L=L1+L2-1;
n1=zeros(L,2*L3);n2=zeros(L,2*L3);n3=zeros(L,2*L3);n4=zeros(L,2*L3);%n1代表A1*B1，n2代表A1*B2，n3代表A2*B1，n4代表A2*B2
  %% 第一部分
  for i=1:L1
      for j=1:L2
          for k=1:L3
              n1(i+j-1,2*k-1)=n1(i+j-1,2*k-1)+1/2*A(i,2*k-1)*B(j,2*k-1);%第一个式子代表积化和差的第一项
              n1(abs(i-j)+1,2*k-1)=n1(abs(i-j)+1,2*k-1)+1*A(i,2*k-1)*B(j,2*k-1)/2;%第二个式子代表积化和差的第二项
          end
      end
  end
  %% 第二部分
  for i=1:L1
      for j=2:L2   %为什么这里是从2开始，因为sin0肯定等于0，取sin0的那一个矩阵的第一行不产生系数
          for k=1:L3
              n2(i+j-1,2*k)=n2(i+j-1,2*k)+1/2*A(i,2*k-1)*B(j,2*k);
              n2(abs(i-j)+1,2*k)=n2(abs(i-j)+1,2*k)-extractsym(i-j)*A(i,2*k-1)*B(j,2*k)/2;
          end
      end
  end
  %% 第三部分
  for i=2:L1
      for j=1:L2
          for k=1:L3
              n3(i+j-1,2*k)=n3(i+j-1,2*k)+1/2*A(i,2*k)*B(j,2*k-1);
              n3(abs(i-j)+1,2*k)=n3(abs(i-j)+1,2*k)+extractsym(i-j)*A(i,2*k)*B(j,2*k-1)/2;
          end
      end
  end
  %% 第四部分
  for i=2:L1
      for j=2:L2
          for k=1:L3
              n4(i+j-1,2*k-1)=n4(i+j-1,2*k-1)-1/2*A(i,2*k)*B(j,2*k);
              n4(abs(i-j)+1,2*k-1)=n4(abs(i-j)+1,2*k-1)+1*A(i,2*k)*B(j,2*k)/2;
          end
      end
  end
%实际上不需要四个部分，一个循环就可以求出，但是这样编程便于修改和纠正
y=n1+n2+n3+n4;

%% 以下是可能会用得到的代码
% A=[1 2 3 4 5 6 7 8;2 2 3 4 5 6 7 8 ;3 2 3 4 5 6 7 8;4 2 3 4 5 6 7 8]
% B=[1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1]
% y=matrixmul_2(A,B)
% A=[0 0;1 1]
% B=[1 0;0 0;0 1]

