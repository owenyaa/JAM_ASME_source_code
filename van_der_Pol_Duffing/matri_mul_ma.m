function y=matri_mul_ma(A,B)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %�����Ǻ���ϵ������ͷ�����˵ĺ���A
  %A�Ƿ���N_dof*N_dof),B������ϵ������n*2N_dof)
  %�����������ϵ������n*2N_dof)
  %% �ú���ʵ�ֵķ�������
  %
%% ʵ��
[n,N_dof]=size(B);
N_dof=N_dof/2;
for i=1:N_dof
    X(i).part=B(:,2*i-1:2*i);
end
for i=1:N_dof
    solution(i).part=A(i,1)*X(1).part;
end
for i=1:N_dof
    for j=1:N_dof-1
        solution(i).part=matrixsum_2(solution(i).part,A(i,j+1)*X(j+1).part);
    end
end
y=zeros(n,2*N_dof);
for i=1:N_dof
    y(:,2*i-1:2*i)=solution(i).part;
end