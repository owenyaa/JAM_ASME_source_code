function y=matrixsum_2(A,B)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %�������Ǻ���ϵ������cosnt,sinntǰ���ϵ�����͵ĺ���
  %A��B����Ҫ���ϵ������Ҳ��n*2�ľ���
  %% �ú���ʵ�ֵķ�������
  %���ж�A��B������ĳ���
  %���½�һ�����󣬳���ȡA,B���ȵ����ֵ
  %ȡ���ȴ�ľ���ΪA��С��ΪB
  %�µľ�����ͷ��ÿһ��Ԫ�ض�����AB��Ԫ��֮��

%% ʵ��
  %�������ά��
  L1=size(A,1);L2=size(B,1);
  L3=size(A,2)/2;
  L=max(L1,L2);
  solution=zeros(L,2*L3);  
  
  %��������
  if L1 < L2
      A1=B;
      B=A;
      A=A1;
      L2=L1;
      L1=L;
  end 
  
  %%�������Ԫ��
  for i=1:L2
      for j=1:2*L3
          solution(i,j)=A(i,j)+B(i,j);
      end
  end
  for i=L2+1:L1
      for j=1:2*L3
          solution(i,j)=A(i,j);
      end
  end
 for i=1:L3
 solution(1,2*i)=0;
 end
  y=solution;
  %% �����ǿ��ܻ��õõ��Ĵ���
  %A=[1 2 3 4 5 6 7 8;2 2 3 4 5 6 7 8 ;3 2 3 4 5 6 7 8;4 2 3 4 5 6 7 8]
  %B=[1 1 1 1 1 1 1 1;1 1 1 1 1 1 1 1]
  %y=matrixsum_2(A,B)