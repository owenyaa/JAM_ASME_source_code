function y=arrange_column_inv(A)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %�ú�����һ������ת��Ϊ���о��󣨲��������Ǿ�����ʽ��
  %������������������о���
  %% �ú���ʵ�ֵķ�������
  %�������ɶ�ƽ��������
  global N_dof
  L1=size(A,1);
  n=L1/N_dof;
  solution=zeros(n,N_dof);
  for i=1:N_dof
      solution(1:end,i)=A(1+(i-1)*n:i*n,:);
  end
  y=solution;