function y=matrixmul(every)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %���������Ǻ���ϵ������cosnt,sinntǰ���ϵ������˵ĺ���
  %% �ú���ʵ�ֵķ�������
  %every�е�Ԫ��һ��һ�������
%% ʵ��
  N=size(every,2);
  S=every(1).matrix;
  for i=2:N
      S=matrixmul_2(S,every(i).matrix);
  end
  y=S;