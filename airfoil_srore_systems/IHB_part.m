function part=IHB_part(parameter_a,w0)
%% ˵��
  %% �ú������ܵļ򵥽��ܣ�
  %ֻ����������Լ���
  %�����ģ�ͲвJa��Jw�����Ǿ���
  %M,C,K����������նȾ���F�Ƿ�������(Ҫ�Լ���ǰ�������������
  %����в�����������Ja�Ǿ���Jw������������arrange_column�����ֵ
  %% �ú���ʵ�ֵķ�������
  %�����еĲ������뷽������⣬Ȼ�����
%% ʵ��
global Q N_harm 
[M,C,K,N_dof]=MCK(Q);
da=zeros(N_harm,2*N_dof);
a_1=matrixder(parameter_a,1);a_2=matrixder(parameter_a,2);
F=nonterm(w0,parameter_a,da);
%% �в�
R(1).matrix=matri_mul_ma(M,-w0^2*a_2); R(2).matrix=matri_mul_ma(C,-w0*a_1); R(3).matrix=matri_mul_ma(K,-parameter_a); R(4).matrix=-F(1).part;
R_matrix=matrixsum(R);
R_vector=matrixtra(R_matrix(1:N_harm,:));
R_vector=arrange_column(R_vector);
%% Ja
Ja=jacobi(parameter_a,w0);
%% Jw
Jw(1).matrix=matri_mul_ma(M,2*w0*a_2); Jw(2).matrix=matri_mul_ma(C,a_1); Jw(3).matrix=F(3).part;
Jw=matrixsum(Jw); 
Jw=matrixtra(Jw(1:N_harm,:));
Jw=arrange_column(Jw);
%% ���
part(1).vector=R_vector; part(2).vector=Ja; part(3).vector=Jw;