function y=matrixsum(every)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %�������Ǻ���ϵ������cosnt,sinntǰ���ϵ�����͵ĺ���
  %every����Ҫ���ϵ���������飬��ͷ�ж�����󣬶���n*2�ľ���
  %% �ú���ʵ�ֵķ�������
  %���ж�every�о���ĳ���
  %���½�һ�����󣬳���ȡ���ȵ����ֵ
  %��������ϵ��������ӵ�С����һֱ���

%% ʵ��
  %�������ά��
  N=size(every,2);
  len=zeros(N,1);
  for i=1:N
      len(i)=size(every(i).matrix,1);
      L1=max(len);
  end
  L2=size(every(1).matrix,2);
  S=zeros(L1,L2);
  %�������Ӻ������
  for i=1:N
      S=matrixsum_2(S,every(i).matrix);
  end
  
  y=S;