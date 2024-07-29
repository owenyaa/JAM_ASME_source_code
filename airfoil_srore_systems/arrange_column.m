function y=arrange_column(A)
%% 说明
  %% 该函数功能的简单介绍：
  %该函数将多列矩阵排为一列
  %% 该函数实现的方法介绍
  %第N+1列放在N列下面
  [L1,L2]=size(A);
  S=zeros(L1*L2,1);
  for i=1:L2
      S(L1*(i-1)+1:L1*(i),1)=A(:,i);
  end
  y=S;