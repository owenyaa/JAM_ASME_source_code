function y=matri_sum_con(A,B)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %�����Ǻ���ϵ������ͷ�����˵ĺ���A
  %A�ǳ���,B������ϵ������n*2N_dof)
  %�����������ϵ������n*2N_dof)
  %% �ú���ʵ�ֵķ�������
  %
%% ʵ��
B(1,1)=B(1,1)+A;
y=B;