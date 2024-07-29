%% 说明
  %% 该程序功能的简单介绍：
  %用IHB法求解非线性方程组（自激问题）
  %% 该程序实现的方法介绍
  %给定初值，通过循环求解增量的方法迭代求解方程的自激运动及其频率
  %具体理论方法参见IHB法，这里面只给出计算机实现过程中的重点
%% 实现
clear
clc
syms t a01 a11 b11 a21 b21 a02 a12 b12 a22 b22 a03 a13 b13 a23 b23
a01
a11
b11
a21
b21
a02
a12
b12
a22
b22 
a03 
a13 
b13 
a23
b23

x1=a01+a11*cos(1*t)+a21*cos*(2*t)+b11*sin(1*t)+b21*sin(2*t);dx1_1=diif(x1,1);dx1_2=diff(x1,2);
x2=a02+a11*cos(1*t)+a12*cos*(2*t)+b12*sin(1*t)+b22*sin(2*t);dx2_1=diif(x2,1);dx2_2=diff(x2,2);
x3=a03+a13*cos(1*t)+a13*cos*(2*t)+b13*sin(1*t)+b23*sin(2*t);dx3_1=diif(x3,1);dx3_2=diff(x3,2);


