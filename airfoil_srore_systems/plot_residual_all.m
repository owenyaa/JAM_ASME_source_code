clear;close all;
load TIHB_N20_LCO1.mat;
dt=0.01;
Tdata=0:dt:1000;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm-1
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(i)*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*(i)*Harm_parameter_a(i,2*j)*cos((i)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(i))^2*Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)-(w0*(i))^2*Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
    end
    x(j,:)=x(j,:)+parameter_a(2*j-1,1);
end

k_h3=0;k_alpha3=10;k_beta3=0;
F=[k_h3*x(1,:).^3;k_alpha3*x(2,:).^3;k_beta3*x(3,:).^3];
residual=M*ddx+C*dx+K*x+F;

subplot(3,1,1);
plot(Tdata,residual(2,:),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
clear;
load IHB_N20_LCO1.mat;
dt=0.01;
Tdata=0:dt:1000;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm-1
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(i)*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*(i)*Harm_parameter_a(i,2*j)*cos((i)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(i))^2*Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)-(w0*(i))^2*Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
    end
    x(j,:)=x(j,:)+parameter_a(2*j-1,1);
end

k_h3=0;k_alpha3=10;k_beta3=0;
F=[k_h3*x(1,:).^3;k_alpha3*x(2,:).^3;k_beta3*x(3,:).^3];
residual=M*ddx+C*dx+K*x+F;

hold on;
plot(Tdata,residual(2,:),'r--','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

clear;
load TIHB_N20_LCO2.mat;
dt=0.01;
Tdata=0:dt:1000;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm-1
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(i)*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*(i)*Harm_parameter_a(i,2*j)*cos((i)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(i))^2*Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)-(w0*(i))^2*Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
    end
    x(j,:)=x(j,:)+parameter_a(2*j-1,1);
end

k_h3=0;k_alpha3=10;k_beta3=0;
F=[k_h3*x(1,:).^3;k_alpha3*x(2,:).^3;k_beta3*x(3,:).^3];
residual=M*ddx+C*dx+K*x+F;

subplot(3,1,2);
plot(Tdata,residual(2,:),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
clear;
load IHB_N20_LCO2.mat;
dt=0.01;
Tdata=0:dt:1000;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm-1
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(i)*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*(i)*Harm_parameter_a(i,2*j)*cos((i)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(i))^2*Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)-(w0*(i))^2*Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
    end
    x(j,:)=x(j,:)+parameter_a(2*j-1,1);
end

k_h3=0;k_alpha3=10;k_beta3=0;
F=[k_h3*x(1,:).^3;k_alpha3*x(2,:).^3;k_beta3*x(3,:).^3];
residual=M*ddx+C*dx+K*x+F;

hold on;
plot(Tdata,residual(2,:),'r--','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);

clear;
load TIHB_N20_LCO3.mat;
dt=0.01;
Tdata=0:dt:1000;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm-1
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(i)*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*(i)*Harm_parameter_a(i,2*j)*cos((i)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(i))^2*Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)-(w0*(i))^2*Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
    end
    x(j,:)=x(j,:)+parameter_a(2*j-1,1);
end

k_h3=0;k_alpha3=10;k_beta3=0;
F=[k_h3*x(1,:).^3;k_alpha3*x(2,:).^3;k_beta3*x(3,:).^3];
residual=M*ddx+C*dx+K*x+F;

subplot(3,1,3);
plot(Tdata,residual(2,:),'k-','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);
clear;
load IHB_N20_LCO3.mat;
dt=0.01;
Tdata=0:dt:1000;
x=zeros(N_dof,length(Tdata));dx=zeros(N_dof,length(Tdata));ddx=zeros(N_dof,length(Tdata));
Harm_parameter_a=parameter_a(2:end,:);
for j=1:N_dof
    for i=1:N_harm-1
        x(j,:)=x(j,:)+Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)+Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
        dx(j,:)=dx(j,:)-w0*(i)*Harm_parameter_a(i,2*j-1)*sin(i*w0*Tdata)+w0*(i)*Harm_parameter_a(i,2*j)*cos((i)*w0*Tdata);
        ddx(j,:)=ddx(j,:)-(w0*(i))^2*Harm_parameter_a(i,2*j-1)*cos((i)*w0*Tdata)-(w0*(i))^2*Harm_parameter_a(i,2*j)*sin((i)*w0*Tdata);
    end
    x(j,:)=x(j,:)+parameter_a(2*j-1,1);
end

k_h3=0;k_alpha3=10;k_beta3=0;
F=[k_h3*x(1,:).^3;k_alpha3*x(2,:).^3;k_beta3*x(3,:).^3];
residual=M*ddx+C*dx+K*x+F;

hold on;
plot(Tdata,residual(2,:),'r--','LineWidth',1);
set(gca,'FontName','Times New Roman','FontSize',15,'LineWidth',1.5);


