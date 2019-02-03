////////////////////////////////////////////////////////////////////////////////////
//�ύX��CEE(Christiano, Eichenbaum and Evans, 2005)���f��
//���K�̓}�N��DSGE���f��(�J���u���[�V����)
//2018�N12��7��

//Please note that the following code is only for "Dynare for MATLAB".
//You can use/rewrite this code without my permission.
//I will not be responsible for any damage/liability induced by running this code.
/////////////////////////////////////////////////////////////////////////////////////

//�p�[�g�P�F�ϐ��̒�`///////////////////////////////////

//�����ϐ�
var y c co cr lambda i q k kg n w rk b r pi g gi tau z zeta ucc uii unn;

//�O���ϐ�
varexo eg ez et em eq ec ei egi ew en;

//�p�����[�^
parameters beta theta varphi h omega alpha nu deltap deltag lw rw zw;
parameters kappa eta psi phi phih phir phipi phiy phit phiyt phibt phiig;
parameters rhog rhoz rhoc rhoi rhogi rhon bybar gybar;

//�p�[�g�Q�F�J���u���[�g����p�����[�^�̒�`(�]��2011)//////////////////////////////////

h    = 0.7742 ; %����҂̏K���`���̓x����
beta = 0.996; %������
theta = 1.5; %����̑�֒e�͐�
varphi = 2; %�J�������̑�֒e�͐�
alpha = 0.33; %���{���z��
omega   =     0.4230   ;   %�񃊃J�[�h�ƌv�̊���
deltap = 0.06; %���Ԏ��{���՗�
deltag = 0.04; %�Љ�{���՗�
eta = 0.75; %Calvo�p�����[�^�[
kappa = 7; %�����̒����R�X�g
psi = 11; %�}�[�N�A�b�v��
zw  = 0.3760      ;  %�������œK���ł��Ȃ��m��
lw  = 0.2; %�����}�[�N�A�b�v�V���b�N
rw  = 0.2; %�����̉ߋ��̃C���t�����ɑ΂���A���x����
nu      =     0.2520    ;   %�Љ�{�̐��Y�͌���
phi      =     0.6509    ;   %���Z����V���b�N�̎�����
phih    = 0.1231  ; %�����̉ߋ��̃C���t�����ɑ΂���A���x����(�n�C�u���b�hNKPC)
phipi     =    1.7423     ;  %�e�C���[���[���̃C���t�����̍�
phir      =   0.0835       ; %�e�C���[���[���̖��ڋ����̎�����
phiy      =   0.3272       ; %�e�C���[���[����GDP��
phit      =   0.6382       ; %�ېŃ��[���̑O����
phiig      = 0.0050        ;  %���������V���b�N�̑O��GDP��
phibt     =   0.1794       ; %�ېŃ��[���̍��̍�
phiyt      =  0.0070      ; %�ېŃ��[����GDP�̍�
rhog   =      0.5523      ; %��������V���b�N�̎�����
rhogi        =0.4166      ; %���������V���b�N�̎�����
rhoz    =     0.8765       ; %�Z�p�V���b�N�̎�����
rhoc         =0.4737      ; %�I�D�V���b�N�̎�����
rhoi         =0.49      ; %�����V���b�N�̎�����
rhon         =0.51       ; %�J�������V���b�N�̎�����

bybar = 2.35; %���ɂ����鍑�̑�GDP��
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
q=pi(+1)-r+((1-deltap)/(1-deltap+rkbar))*q(+1)+(rkbar/(1-deltap
+rkbar))*rk(+1)+eq; %�g�[�r����q
//pi=beta*pi(+1)+(((1-eta)*(1-beta*eta))/(eta))*((1-alpha)*w+alpha*rk-
//nu*kg(-1)-z); %�W��NKPC
//�n�C�u���b�hNKPC///////////////////////////////////////
pi=(beta/(1+beta*phih))*pi(+1)+(phih/
1+beta*phih)*pi(-1)+(1/1+beta*phih)*(((1-eta)*(1-beta*eta))/(eta))*((1-
alpha)*w+alpha*rk-nu*kg(-1)-z); %�n�C�u���b�hNKPC
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
unn=rhon*unn(-1)+en; %�J�������V���b�N
end;

//�p�[�g�S�F�V���b�N�A�V�~�����[�V����/////////////////////////////////////
shocks;
var et = 1; %���ŃV���b�N
var em = 1; %���Z�ɘa�V���b�N
var ez = 1; %�Z�p�V���b�N
var eg = 1; %���{����V���b�N
var eq = 1; %�،��v���~�A���V���b�N
var ec = 1; %�I�D�V���b�N
var ei = 1; %�����V���b�N
var ew = 1; %�����}�[�N�A�b�v�V���b�N
var en = 1; %�J�������V���b�N
var egi = 1; %���������V���b�N
end;

check; resid(1); steady; %BK�����̊m�F�A����Ԓl 

stoch_simul(irf=40) y c i n k b w r pi ; %�m���I�V�~�����[�V����(40��)
//close all; %//�������Ǝ����I�ɃO���t��������悤�ɂȂ�(�v�Z���������)�B