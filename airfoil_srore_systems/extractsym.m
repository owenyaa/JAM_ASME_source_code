function y=extractsym(A)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %���һ�����ķ��ţ��������Ǹ���
  %A��Ҫ��ȡ���ŵ���ֵ
  %% �ú���ʵ�ֵķ�������
  %����������Լ��ľ���ֵ������Ƿ�0���Ļ�
%% ʵ��
if A==0
    %'The value entered is 0, please check'
    y=0;
else
    y=A/abs(A);
end
    