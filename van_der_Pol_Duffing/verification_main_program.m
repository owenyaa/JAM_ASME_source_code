%% ˵��
  %% �ó����ܵļ򵥽��ܣ�
  %��IHB���������Է����飨���⼤�������������⣩
  %% �ó���ʵ�ֵķ�������
  %������ֵ��ͨ��ѭ����������ķ���������ⷽ�̵��˶�����Ƶ��
  %�������۷����μ�IHB����������ֻ���������ʵ�ֹ����е��ص�
%% ʵ��
clear
clc
close all
tic
global w0 a n M C K N_dof
%w0���Լ��񶯵�Ƶ��Ҳ��δ֪����a�����Ǽ�����ϵ��Ҳ��δ֪����n�����õ�г������,Q�Ƿ��٣���һ���õ�����N_dof�����ɶȸ���
n=20;iter=0;
[M,C,K,N_dof]=MCK;
w0=0.3;
% [a,~]=get_initial(n);a(5:end,:)=zeros(n-5+1,2*N_dof);

a=zeros(n,2*N_dof);
a(1,:)=[0,0];
a1=0.3;a2=-0.9;
a(2,:)=[a1,a2];

%% ѭ������������Ͻӽ���ȷ��
er=1;tol=1e-6;
dw=0;
while(er>=tol&&iter<1001)
    part=IHB_part(M,C,K);
    R=part(1).vector; Ja=part(2).vector; Jw=part(3).vector;
    er=norm(R)
    da=Ja\R;
    %     [u,s,v]=csvd(Ja);
    %     lambda=l_curve(u,s,R);
    %     da=tikhonov(u,s,v,R,lambda);
    da=arrange_column_inv(da);
    da=matrixtra_inv(da);
    w0=w0+dw;
    a=a+da;
    iter=iter+1;
end
toc