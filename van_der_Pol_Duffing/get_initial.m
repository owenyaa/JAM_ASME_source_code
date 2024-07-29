%% 说明
%% 该程序功能的简单介绍：
%用最小二乘原理获得IHB法的初值
%% 该程序实现的方法介绍
%通过ODE45来求解该系统，将获得的数值解进行最小二乘，可以求得谐波系数
%% 实现
function [a,w]=get_initial(n,w0)
global Q
% n=20;
% Q=6.2851*0.75;
[~,~,~,N_dof]=MCK;
a=zeros(n,2*N_dof);
Tend=10000;
dt=0.1;

Tspan=0:dt:Tend;
% y0=ones(2*N_dof,1);
y0=[1,1];
options=odeset('RelTol',1e-10);
[T,X]=ode45(@(t,x) odefun(w0,t,x), Tspan, y0, options);

% figure;
% plot(X(end-500/dt:end,1),X(end-500/dt:end,1+N_dof),'r.','LineWidth',1)
% hold on
% plot(X(end-500/dt:end,2),X(end-500/dt:end,2+N_dof),'r.','LineWidth',1)
% hold on
% plot(X(end-500/dt:end,3),X(end-500/dt:end,3+N_dof),'r.','LineWidth',1)



% [value,locate]=findpeaks(X(end-500/dt:end,1));
% period=(locate(end)-locate(end-1))*dt;
% w0=2*pi/period
w=w0;
period=2*pi/w0;
n_dot=floor(2*period/dt);

% dt=Tspan(2)-Tspan(1);
% data=X(end-2000/dt:end,1);
% N=length(data);
% N_fft=2^16;
% Y=fft(data,N_fft);
% Pyy=Y.*conj(Y)/(N_fft*N);
% f=2*pi/dt*(0:N_fft/2)/N_fft;
% figure
% plot(f,(Pyy(1:(N_fft/2+1))))
% title('Frequency content of y')
% xlabel('frequency (Hz)')

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

for i=1:N_dof
    Y(i).part=X(end-n_dot:end,i);
end

t=T(end-n_dot:end,1);

for j=1:N_dof
    CS(j).part=zeros(size(t,1),2*n-1);
    CS(j).part(:,1)=1;
end
for j=1:N_dof
    for i=2:n
        CS(j).part(:,i)=cos((i-1)*t*w0);
        CS(j).part(:,i+n-1)=sin((i-1)*t*w0);
    end
end

for i=1:N_dof
    A(i).part=(transpose(CS(i).part)*CS(i).part)\(transpose(CS(i).part)*Y(i).part);
    B(i).part=matrixtra_inv(A(i).part);
end
for i=1:N_dof
    a(:,2*i-1:2*i)=B(i).part;
end
end


function y=odefun(w0,t,x)

[M,C,K,N_dof]=MCK;
beta=0.3;f=1;r=0.2;

F=beta*x(1)^3+r*x(1)^2*x(2);
f_ex=f*cos(w0*t);
y=zeros(2*N_dof, 1);
y(1:N_dof)=x(N_dof+1:end);
y(N_dof+1:end)=M\(f_ex-K*x(1:N_dof)-C*x(N_dof+1:end)-F);

end