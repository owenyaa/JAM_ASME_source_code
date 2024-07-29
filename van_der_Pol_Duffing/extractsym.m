function y=extractsym(A)
%% 说明
  %% 该函数功能的简单介绍：
  %获得一个数的符号（是正还是负）
  %A是要提取符号的数值
  %% 该函数实现的方法介绍
  %这个数除以自己的绝对值，如果是非0数的话
%% 实现
if A==0
    %'The value entered is 0, please check'
    y=0;
else
    y=A/abs(A);
end
    