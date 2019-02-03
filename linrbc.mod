////////////////////////////////////////////////////////////////////////////////////
	//���`RBC���f��(�J���u���[�V����)
	//Boxcat(������w�o�ϊw���S�N) 
	//2018�N12��8��
	
	//Please note that the following code is only for "Dynare for MATLAB".
	//You can use/rewrite this code without my permission.
	//I will not be responsible for any damage/liability induced by running this code.
	/////////////////////////////////////////////////////////////////////////////////////
//RBC���f����Dynare�R�[�h
//�p�[�g�P�@�ϐ��̒�`
var c n y r w k i z;%�����ϐ�
varexo e;%�O���ϐ�

//�p�[�g�Q�@�p�����[�^�̐ݒ�
parameters alpha beta theta phi delta rho;%�p�����[�^���`
//�p�����[�^�[�ɒl����
alpha=0.33;%���{���z��
beta=0.996;%������
theta=1.5;%����̑�֒e�͐��̋t��
phi=2;%�J�������̑�֒e�͐��̋t��
delta=0.04;%���{���՗�
rho=0.9;%�Z�p�i���̎�����

//�p�[�g3���f���̋L�q
model(linear);

%����Ԓl�ŋ��߂����l
#rbar=(1+beta*delta-beta)/theta;%r/(1+r-delta);
#iybar=(delta*alpha*beta)/(1+beta*delta-beta);
#cybar=1-iybar;

y=z+alpha*k(-1)+(1-alpha)*n;%���Y�֐�
r=z+(alpha-1)*k(-1)+(1-alpha)*n;%���{�̃����^����
w=z+alpha*k(-1)+alpha*n;%�J�����v�̍œK��
c=c(+1)-rbar*r(+1);%����̃I�C���[������
n=(1/phi)*w-(theta/phi)*c;%�œK�J����������
y=cybar*c+iybar*i;%���s��̋ύt����
k=(1-delta)*k(-1)+delta*i;%���{�̐��ڎ�
z=rho*z(-1)+e;%�Z�p�V���b�N
end;

//�p�[�g4�@���f���̃`�F�b�N
check;
resid(1);
steady;

//�p�[�g5�@�V���b�N
shocks;
var e = 0.01;
end;

//�p�[�g6�@�V�~�����[�V����(100��)
stoch_simul(irf=100)y c i n w k r z;