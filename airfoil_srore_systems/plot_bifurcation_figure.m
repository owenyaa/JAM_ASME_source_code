clear;%close all;

clear;
load new_P1.mat;
% QQ=2.3:0.05:3.2;%最下方线
QQ=4.5:-0.01:2.58;%P1
% QQ=3:-0.05:2.6;%P2
for ii=1:1:length(QQ)
    w01(ii)=every(ii).w;Q1(ii)=every(ii).Q;
    %parameter_a1(ii)=every(ii).parameter_a;
    A_211(ii)=every(ii).A_21;
end
hold on;
plot(Q1,A_211,'r-')

clear;
load new_P2.mat;
QQ=2.59:0.01:3.21;%P2_1
for ii=1:1:length(QQ)
    w01(ii)=every(ii).w;Q1(ii)=every(ii).Q;
%     parameter_a1(ii)=every(ii).parameter_a;
    A_211(ii)=every(ii).A_21;
end
hold on;
plot(Q1,A_211,'b-')

load new_P3.mat;
QQ=2.3:0.01:3.21;%P3
% QQ=4.5:-0.05:2.6;%最上方线
% QQ=3:-0.05:2.6;%P2
for ii=1:1:length(QQ)
    w01(ii)=every(ii).w;Q1(ii)=every(ii).Q;
    %parameter_a1(ii)=every(ii).parameter_a;
    A_211(ii)=every(ii).A_21;
end
figure;
plot(Q1,A_211,'k-')




