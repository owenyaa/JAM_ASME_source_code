function y=TRIseries(n)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %�������n�׵����Ǿ���ĺ���
  %�������������
  %% �ú���ʵ�ֵķ�������
  %��a0,0;cos(t),sin(t);cos(2*t);sin(2*t).....��
%% ʵ��
syms t
S=sym(zeros(n,2));
for i=1:n
    if i==1
        S(i,1)=1;
    else
        S(i,1)=cos((i-1)*t);
    end
end
for i=1:n
    if i==1
        S(i,2)=0;
    else
        S(i,2)=sin((i-1)*t);
    end
end
y=sym(zeros(2*n-1,1));
y(1:n,1)=S(:,1);
y(n+1:end,1)=S(2:end,2);

