//////////////////////////////////////////////////////////////////////////////
//�ύX��CEE(Christiano, Eichenbaum and Evans, 2005)���f��
//���K�̓}�N���v��DSGE���f��
//2018�N12��7��

//Please note that the following code is only for "Dynare for MATLAB".
//You can use/rewrite this code without my permission.
//I will not be responsible for any damage/liability induced by running this code.
////////////////////////////////////////////////////////////////////////////

//�p�[�g�P�F�ϐ��̒�`///////////////////////////////////

//�����ϐ�
var y c co cr lambda i q n k kg w rk b r pi g gi tau z zeta ucc uii unn;
var y_obs c_obs i_obs g_obs pi_obs r_obs w_obs n_obs b_obs;

//�O���ϐ�
varexo ec ei ew en eg egi ez et em eq uy uc ui ug upi ur uw un ub;

//�p�����[�^
parameters beta theta varphi h omega alpha nu deltap deltag lw rw zw;
parameters kappa eta psi phi phih phir phipi phiy phit phiyt phibt phiig;
parameters rhog rhoz rhoc rhoi rhogi rhon bybar gybar;

//�p�[�g�Q�F�J���u���[�g����p�����[�^�̒�`(�]��2011)//////////////////////////////////

beta = 0.996; %������
theta = 1.5; %����̑�֒e�͐�
varphi = 2; %�J�������̑�֒e�͐�
alpha = 0.33; %���{���z��
deltap = 0.06; %���Ԏ��{���՗�
deltag = 0.04; %�Љ�{���՗�
lw=0.2; %�����}�[�N�A�b�v�V���b�N
rw=0.2; %�����̉ߋ��̃C���t�����ɑ΂���A���x����
eta = 0.75; %Calvo�p�����[�^�[
kappa = 7; %�����̒����R�X�g
psi = 11; %�}�[�N�A�b�v��

bybar = 4; %���ɂ����鍑�̑�GDP��
gybar = 0.08; %���ɂ�������������̑�GDP��

//�p�[�g�R�F���f��///////////////////////////////////////

model(linear);

//����Ԓl�̒�`

# rbar = 1/beta; %���e���q���̒���Ԓl
# rkbar=rbar+deltap-1; %���{�̃����^����
# wbar=(((psi-1)/psi)*((((1-alpha)^(1-alpha))*(alpha^(alpha)))/
(rkbar^(alpha))))^(1/(1-alpha)); %���ɂ����������
# nybar=((((1-alpha)*rkbar)/(alpha*wbar)))^(alpha); %���ɂ�����J���ʂ̑�GDP��
# kybar=((((1-alpha)*rkbar)/(alpha*wbar)))^(alpha-1); %���ɂ����鎑�{�̑�GDP��
# iybar=deltap*kybar; %���ɂ��������̑�GDP��
# cybar=1-iybar-gybar; %���ɂ����铊���̑�GDP��

//���^���f��

(1-h)*(1-beta*h)*lambda = -theta*(co-h*co(-1))+(1-h)*ucc
+beta*h*(theta*(co(+1)-h*co)-(1-h)*ucc(+1)); %���J�[�h�I�ƌv�̏���
lambda = lambda(+1)+r-pi(+1);
cr=((wbar*nybar)/cybar)*(w+n)-(1/cybar)*tau; %�񃊃J�[�h�I�ƌv�̏���
c=omega*cr+(1-omega)*co; %�o�ϑS�̂̏���
//��������////////////////////////////////////////////////////
w=(beta/(1+beta))*w(+1)+(1/(1+beta))*w(-1)+(beta/(1+beta))*pi(+1)-
((1+beta*rw)/1+beta)*pi+
(rw/(1+beta))*pi(-1)-(1/(1+beta))*((lw*(1-beta*zw)*(1-zw))/(lw
+(1+lw)*varphi)*zw)*(w-varphi*n+lambda-unn)-ew;
/////////////////////////////////////////////////////////////
i=(1/(1+beta))*i(-1)+(beta/(1+beta))*i(+1)+(kappa/(1+beta))*q+(beta/
(1+beta))*(uii(+1)-uii); %�����֐�
q=pi(+1)-r+((1-deltap)/(1-deltap+rkbar))*q(+1)+(rkbar/(1-deltap+rkbar))*rk(+1)
+eq; %�g�[�r����q
//pi=beta*pi(+1)+(((1-eta)*(1-beta*eta))/(eta))*((1-alpha)*w+alpha*rk-
//nu*kg(-1)-z); %�W��NKPC
//�n�C�u���b�hNKPC///////////////////////////////////////
pi=(beta/(1+beta*phih))*pi(+1)+(phih/1+beta*phih)*pi(-1)+(1/1+beta*phih)*(((1-
eta)*(1-beta*eta))/(eta))*((1-alpha)*w+alpha*rk-nu*kg(-1)-z); %�n�C�u���b�hNKPC
////////////////////////////////////////////////////////
k=(1-deltap)*k(-1)+deltap*i; %���{�̐��ڎ�
y=z+alpha*k(-1)+(1-alpha)*n+nu*kg(-1); %���Y�֐�
y=cybar*c+iybar*i+g +gi; %���s��̋ύt����
r=phir*r(-1)+(1-phir)*(phipi*((pi+pi(-1)+pi(-2)+pi(-3))/4)+phiy*y)+zeta; 
%�e�C���[���[��(�����X���[�W���O)
n-k(-1)=rk-w; %��p�ŏ�������
b=rbar*b(-1)+bybar*rbar*r(-1)-rbar*bybar*pi+gi+g-tau; %���{�̗\�Z����
kg=(1-deltag)*kg(-1)+deltag*(1/gybar)*g; %�Љ�{�̒~�ώ�
tau=phit*tau(-1)+(1-phit)*(phiyt*y+phibt*b(-1))-et; %�������̂���Ő����[��
gi=rhogi*gi(-1)+phiig*y(-1)+egi; %���������V���b�N
g=rhog*g(-1)+eg; %���{����V���b�N
z=rhoz*z(-1)+ez; %�Z�p�V���b�N
zeta=phi*zeta(-1)-em; %���Z����V���b�N
ucc=rhoc*ucc(-1)+ec; %�I�D�V���b�N
uii=rhoi*uii(-1)+ei; %�����V���b�N
unn=rhon*unn(-1)+en; %�J���V���b�N

//�ϑ�������

y_obs=y+uy; %GDP�̊ϑ�������
c_obs=c+uc; %����̊ϑ�������
i_obs=i+ui; %�����̊ϑ�������
g_obs=g+gi+ug; %���{�x�o�̊ϑ�������
pi_obs=pi+upi; %�C���t�����̊ϑ�������
r_obs=r+ur; %���ڗ��q���̊ϑ�������
w_obs=w+uw; %�������̊ϑ�������
n_obs=n+un; %�J���ʂ̊ϑ�������
b_obs=b+ub; %���̊ϑ�������
end;

//�p�[�g�S�F���v����p�����[�^�̎��O���z(�]��2011�A�A��2012�A���c2009�Ȃ�)///////////

estimated_params;

h, normal_pdf, 0.8, 0.15; %����҂̏K���`���̓x����
omega, beta_pdf, 0.3, 0.1; %�񃊃J�[�h�ƌv�̊���
nu, normal_pdf, 0.25, 0.1; %�Љ�{�̐��Y�͌���
zw, beta_pdf, 0.375, 0.1; %�������œK���ł��Ȃ��m��
phi, normal_pdf, 0.22, 0.1; %���Z����V���b�N�̎�����
phipi, normal_pdf, 1.775, 0.2; %�e�C���[���[���̃C���t�����̍�
phir, normal_pdf, 0.02, 0.1; %�e�C���[���[���̖��ڋ����̎�����
phiy, normal_pdf, 0.1, 0.1; %�e�C���[���[����GDP��
phih, beta_pdf, 0.25, 0.1; 
%�����̉ߋ��̃C���t�����ɑ΂���A���x����(�n�C�u���b�hNKPC)
phit, beta_pdf, 0.8, 0.1; %�ېŃ��[���̑O����
phiig, normal_pdf, 0.1, 0.1; %���������V���b�N�̑O��GDP��
phibt, normal_pdf, 0.1, 0.1; %�ېŃ��[���̍��̍�
phiyt, normal_pdf, 0.12, 0.1; %�ېŃ��[����GDP�̍�
rhog, beta_pdf, 0.6, 0.1; %��������V���b�N�̎�����
rhoz, beta_pdf, 0.9, 0.1; %�Z�p�V���b�N�̎�����
rhoc, beta_pdf, 0.5, 0.2; %�I�D�V���b�N�̎�����
rhoi, beta_pdf, 0.5, 0.2; %�����V���b�N�̎�����
rhogi, beta_pdf, 0.5, 0.2; %���������V���b�N�̎�����
rhon, beta_pdf, 0.5, 0.2; %�J���V���b�N�̎�����
stderr et, inv_gamma_pdf, 0.1, inf; %���ŃV���b�N
stderr eg, inv_gamma_pdf, 0.1, inf; %���{����V���b�N
stderr ez, inv_gamma_pdf, 0.1, inf; %�Z�p�V���b�N
stderr em, inv_gamma_pdf, 0.1, inf; %���Z�ɘa�V���b�N
stderr uy, inv_gamma_pdf, 0.1, inf; %�ȉ��A�e�����ϐ��̊ϑ��덷
stderr uc, inv_gamma_pdf, 0.1, inf;
stderr ui, inv_gamma_pdf, 0.1, inf;
stderr ug, inv_gamma_pdf, 0.1, inf;
stderr upi, inv_gamma_pdf, 0.1, inf;
stderr ur, inv_gamma_pdf, 0.1, inf;
stderr uw, inv_gamma_pdf, 0.1, inf;
stderr un, inv_gamma_pdf, 0.1, inf;
stderr ub, inv_gamma_pdf, 0.1, inf;
end;

//�p�[�g�T�F���v///////////////////////////////////////////////////////////

//�ϑ�����ϐ��̎w��
varobs y_obs c_obs i_obs g_obs pi_obs r_obs w_obs n_obs b_obs;

//���v�R�}���h(�J���}���X���[�W���O����)

//estimation(
//datafile=ch8data, %�f�[�^�t�@�C����(1983�`1998�A���{)
//mode_check, %���[�h�`�F�b�N
//mh_replic=100000, %���s��
//mh_nblocks=2, %�`�F�[���̐�
//mh_drop=0.5, mh_jscale=0.37, %��e���ɉe��
//bayesian_irf, %�x�C�Y���� 
//smoother %�J���}���X���[�W���O
//);

//���v�R�}���h(�J���}���X���[�W���O�Ȃ�)//////////////////////////////
estimation(datafile=ch8data, mode_check, mh_replic=10000, mh_nblocks=2, 
mh_drop=0.5, mh_jscale=0.37, bayesian_irf);
////////////////////////////////////////////////////////////////////////


//�����I���ҋύt���̓��o
stoch_simul;

//�q�X�g���J���������s���R�}���h
shock_decomposition(parameter_set=posterior_mean) y_obs;
