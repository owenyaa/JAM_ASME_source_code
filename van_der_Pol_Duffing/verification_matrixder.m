%% 说明
  %% 该程序功能的简单介绍：
  %目的是验证三角函数系数矩阵求导函数是正确的
  %% 该程序实现的方法介绍
  %通过原三角函数求导法则和函数结果进行对比
%% 实现
clc
clear
syms t c0 c1 c2 c3 c4 c5 s0 s1 s2 s3 s4 s5
c0=1;c1=2;c2=3;c3=4;c4=5;c5=6;
s0=6;s1=5;s2=4;s3=3;s4=2;s5=1;
n=0;
x1=c0+c1*cos(t)+c2*cos(2*t)+c3*cos(3*t)+c4*cos(4*t)+c5*cos(5*t)+s1*sin(t)+s2*sin(2*t)+s3*sin(3*t)+s4*sin(4*t)+s5*sin(5*t);

A=[c0 s0 c0 s0;c1 s1 c1 s1;c2 s2 c2 s2;c3 s3 c3 s3;c4 s4 c4 s4;c5 s5 c5 s5];
L=size(A,1);
xn=diff(x1,n);
%% 将x的n次导写成与A相同的矩阵，便于我们进行判断
B=zeros(L,2);
B(1,1)=int(xn,0,2*pi)/(2*pi);
for i=2:L
    B(i,1)=int(xn*cos((i-1)*t),t,0,2*pi)/pi;
    B(i,2)=int(xn*sin((i-1)*t),t,0,2*pi)/pi;
end
B
y=matrixder(A,n)