clear
clc
close all
%% 数值法画分岔图
w0=0.9;
dw=0.01;
N_dof=1;
record=[];
% for i=1:213
%     w0=w0+dw;
Tend=2000;
dt=0.01;

Tspan=0:dt:Tend;
% y0=ones(2*N_dof,1);
% a1=1;a2=1;
a1=0.011914103330285;a2=0.398153445313372;
% a1=2*rand(1)-1;a2=2*rand(1)-1;
y0=[a1,a2];

options=odeset('RelTol',1e-10);
[T,X]=ode45(@(t,x) odefun(t,x,w0), Tspan, y0, options);
period=2*pi/w0;
n_dot=floor(1.1*period/dt);
% record(i,1)=w0;
% record(i,2)=max(X(end-n_dot:end,1));
% i
figure;
plot(Tspan,X(:,1));

figure;
plot(X(end-500/dt:end,1),X(end-500/dt:end,1+N_dof),'r.','LineWidth',1)
% hold on
% plot(X(end-500/dt:end,2),X(end-500/dt:end,2+N_dof),'r.','LineWidth',1)
% hold on
% plot(X(end-500/dt:end,3),X(end-500/dt:end,3+N_dof),'r.','LineWidth',1)

dt=Tspan(2)-Tspan(1);
data=X(end-2000/dt:end,1);
N=length(data);
N_fft=2^16;
Y=fft(data,N_fft);
Pyy=Y.*conj(Y)/(N_fft*N);
f=2*pi/dt*(0:N_fft/2)/N_fft;
figure
plot(f,(Pyy(1:(N_fft/2+1))))
axis([0 3 -inf inf])
title('Frequency content of y')
xlabel('frequency (Hz)')

% figure;
% plot(X(end-n_dot:5:end,1),X(end-n_dot:5:end,1+N_dof),'r.','LineWidth',1)
% hold on
% plot(X(end-n_dot:5:end,2),X(end-n_dot:5:end,2+N_dof),'b.','LineWidth',1)
% hold on
% plot(X(end-n_dot:5:end,3),X(end-n_dot:5:end,3+N_dof),'k.','LineWidth',1)


% figure;
% plot(T(end-n_dot:5:end,1),X(end-n_dot:5:end,1),'r.','LineWidth',1)
% hold on
% plot(T(end-n_dot:5:end,1),X(end-n_dot:5:end,2),'b.','LineWidth',1)
% hold on
% plot(T(end-n_dot:5:end,1),X(end-n_dot:5:end,3),'k.','LineWidth',1)
% end
% figure
% plot(record(:,1),record(:,2),'K.','MarkerSize',12)
function y=odefun(t,x,w0)

[M,C,K,N_dof]=MCK;
beta=0.3;f=1;r=0.2;

F=beta*x(1)^3+r*x(1)^2*x(2);
f_ex=f*cos(w0*t);
y=zeros(2*N_dof, 1);
y(1:N_dof)=x(N_dof+1:end);
y(N_dof+1:end)=M\(f_ex-K*x(1:N_dof)-C*x(N_dof+1:end)-F);

end