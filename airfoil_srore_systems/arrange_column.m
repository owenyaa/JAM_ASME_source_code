function y=arrange_column(A)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %�ú��������о�����Ϊһ��
  %% �ú���ʵ�ֵķ�������
  %��N+1�з���N������
  [L1,L2]=size(A);
  S=zeros(L1*L2,1);
  for i=1:L2
      S(L1*(i-1)+1:L1*(i),1)=A(:,i);
  end
  y=S;