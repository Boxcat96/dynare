////////////////////////////////////////////////////////////////////////////////////
//�Ő�����RBC���f��(non-linear)
//Boxcat
//2018�N12��11��

//Please note that the following code is only for "Dynare for MATLAB".
//You can use/rewrite this code without my permission.
//I will not be responsible for any damage/liability induced by running this code.
/////////////////////////////////////////////////////////////////////////////////////

//�ϐ��̒�`
var y c i n k rk w z g lambda; %�����ϐ�
varexo tauc tauk taun ez;      %�O���ϐ�

//�p�����[�^�̒�`
parameters beta theta varphi alpha delta rhoz;
parameters yss css iss nss kss rkss wss zss lambdass;

beta = 0.996;   %������
theta=1.5;      %����̑�֒e�͐��̋t��
varphi = 2.0;   %�J�������̑�֒e�͐��̋t��
alpha = 0.33;   %���{���z��
delta = 0.04;   %���{���՗�
rhoz = 0.9;     %�Z�p�V���b�N�̎�����

//���f��(����`)
model;
c^(-theta) = (1+tauc)*lambda;       %����̃I�C���[������
n^varphi = (1-taun)*w*lambda;       %�J���̍œK������
lambda = beta * lambda(+1)*(rk(+1)*(1-tauk)+1-delta); %���O�����W���搔
k = (1-delta)*k(-1) + i;            %���{�̐��ڎ�
y = z*(k(-1)^alpha)*n^(1-alpha);    %���Y�֐�
rk = z*(alpha)*(n/k(-1))^(1-alpha); %���{�̃����^���R�X�g
w = z*(1-alpha)*(k(-1)/n)^(alpha);  %������
y = c + i + g;                      %���s��̋ύt����
g = tauc*c + tauk*rk*k(-1) + taun*w*n;  %�Ŏ�
log(z) = rhoz * log(z(-1)) + ez;    %�Z�p�V���b�N
end;

//����Ԓl�̒�`(ss �� steady state �̈Ӗ��ł�)
css=1;
nss=1;
yss=3;
rkss=0.05;
wss=1;
kss=10;
iss=1;
zss = 1;
lambdass = 1;

//�����l�̐ݒ�
initval;
tauc = 0.05;  %�������_�ł̏����
tauk = 0.2;   %�������_�ł̎��{��
taun = 0.2;   %�������_�ł̏�����
c=css;
n=nss;
y=yss;
w=wss;
k=kss;
i=iss;
rk=rkss;
z=zss;
lambda = lambdass;
ez = 0;
end;
steady;%�����̒���Ԃ̌v�Z

//�I�[�l�̐ݒ�
endval;
tauc = 0.1;   %�I�[���_�ł̏����
tauk = 0.2;   %�I�[���_�ł̎��{��
taun = 0.2;   %�I�[���_�ł̏�����
c=css;
n=nss;
y=yss;
w=wss;
k=kss;
i=iss;
rk=rkss;
z=zss;
lambda = lambdass;
ez = 0;
end;
steady;%�V���b�N��̒���Ԃ̌v�Z

check;%BK�����𖞂������m�F

shocks;
var tauc;         %�V���b�N��^����O���ϐ�
periods 0:5, 5:10;     %�V���b�N�̕ύX���
values 0, 0.08;         %�V���b�N�̑傫��(0�Ȃ�A�i�E���X�����g����)
end;

simul(periods=100); %�V���b�N(����)

//����Ԃ���̘������̌v�Z
ctild=(c-c(1))/c(1);
ntild=(n-n(1))/n(1);
ytild=(y-y(1))/y(1);
rktild=(rk-rk(1))/rk(1);
wtild=(w-w(1))/w(1);
ktild=(k-k(1))/k(1);
itild=(i-i(1))/i(1);
gtild=(g-g(1))/g(1);

//�}�̕`��(����Ԃ���̘������̃C���p���X����)
figure(1)
subplot(2,2,1)
plot(0:100, ytild(1:101)); title('GDP')
subplot(2,2,2)
plot(0:100, ctild(1:101)); title('����')
subplot(2,2,3)
plot(0:100, itild(1:101)); title('����')
subplot(2,2,4)
plot(0:100, wtild(1:101)); title('����')

figure(2)
subplot(2,2,1)
plot(0:100, ntild(1:101)); title('�J��')
subplot(2,2,2)
plot(0:100, rktild(1:101)); title('���{�̃����^���R�X�g')
subplot(2,2,3)
plot(0:100, ktild(1:101)); title('���{�~��')
subplot(2,2,4)
plot(0:100, gtild(1:101)); title('�Ŏ�')

//close all;