function y=matrixder(A,n)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %���������Ǻ���ϵ��������󵼣�cosnt,sinntǰ���ϵ����
  %A����Ҫ�󵼵�ϵ������,n���󵼵Ĵ���
  %���������Ķ��ǣ�г����*2�����ɶȣ�ά�ȵľ���
  %% �ú���ʵ�ֵķ�������
  %���Ǻ����󵼻��ж����ϵ������һ�ε��ͻ���һ��n����cos��-sin��sin��cos
  %ͨ��n/2���������ж�������ε�����ż�ε�
  %ż���е�������ͬ����ڶ�����ͬ���������е�������ͬ�����һ����ͬ��
 %% ʵ��
 %ѭ���������ϵ��
 [L1,L2]=size(A);
 solution=zeros(L1,L2);
 remainder=mod(n,2); 
 for i=1:L1
     for j=1:L2
         if i==1 && n>=1
             solution(i,j)=0;
         end
         if remainder==1
             solution(i,j)=(i-1)^n*(-1)^((n-3+2*(2-mod(j,2)))/2)*A(i,j+(-1)^(mod(j,2)+1));
         else
             solution(i,j)=(i-1)^n*(-1)^(n/2)*A(i,j);
         end
     end
 end
 for i=1:L2/2
 solution(1,2*i)=0;
 end
 y=solution;