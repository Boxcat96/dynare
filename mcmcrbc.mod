////////////////////////////////////////////////////////////////////////////////////
//RBC���f��(�x�C�Y���v) Boxcat
//2018�N12��8��

//Please note that the following code is only for "Dynare for MATLAB".
//You can use/rewrite this code without my permission.
//I will not be responsible for any damage/liability induced by running this code.
/////////////////////////////////////////////////////////////////////////////////////
//�p�[�g�P�@�ϐ��̒�`
var c n y r w k i z y_obs c_obs i_obs w_obs n_obs;%�����ϐ�
varexo e uy uc ui uw un;%�O���ϐ�

//�p�[�g�Q�@�p�����[�^�̐ݒ�
parameters alpha beta theta phi delta rho;%�p�����[�^���`
//�p�����[�^�[�ɒl����
beta=0.996;%������
delta=0.04;%���{���՗�

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

y_obs=y+uy; %GDP�̊ϑ�������
c_obs=c+uc; %����̊ϑ�������
i_obs=i+ui; %�����̊ϑ�������
w_obs=w+uw; %�������̊ϑ�������
n_obs=n+un; %�J���ʂ̊ϑ�������
end;

//�p�[�g4�@��������p�����[�^�̎��O���z
estimated_params;
alpha, beta_pdf, 0.33, 0.1;%���{���z��
theta, gamma_pdf, 1.5, 0.1;%����̑�֒e�͐��̋t��
phi, gamma_pdf, 2, 0.1;%�J�������̑�֒e�͐��̋t��
rho, beta_pdf, 0.85, 0.1;%�Z�p�V���b�N�̎�����
stderr e, inv_gamma_pdf, 0.1, inf;
stderr uy, inv_gamma_pdf, 0.1, inf;
stderr uc, inv_gamma_pdf, 0.1, inf;
stderr ui, inv_gamma_pdf, 0.1, inf;
stderr uw, inv_gamma_pdf, 0.1, inf;
stderr un, inv_gamma_pdf, 0.1, inf;
end;
%�`�F�b�N�R�}���h�͎g���Ă͂Ȃ�Ȃ��B

%�p�[�g5�@�V�~�����[�V����
varobs y_obs c_obs i_obs w_obs n_obs;%�ǂݎ��ϑ��f�[�^

estimation(datafile=rbcdata, mode_check, mh_replic=100000, mh_nblocks=2,
mh_drop=0.5, mh_jscale=0.73, bayesian_irf);

stoch_simul;