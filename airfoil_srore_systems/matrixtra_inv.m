function y=matrixtra_inv(A)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %�����������任Ϊ����ϵ������
  %% �ú���ʵ�ֵķ�������
  %ȡ���������ȵ�һ�루Ҫ����ȡ�������ڶ���Ҫ��һ��0
%% ʵ��

[L1,L2]=size(A);
L1=ceil(L1/2);
solution=zeros(L1,2*L2);
for i=1:L2
    solution(:,2*i-1)=A(1:L1,i);
    solution(2:end,2*i)=A(L1+1:end,i);
end
y=solution;
%% �����ǿ����õĵ��Ĵ���
% y=matrixtra_inv(matrixtra(A))