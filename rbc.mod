////////////////////////////////////////////////////////////////////////////////////
	//����`RBC���f��(�J���u���[�V����)
	//Boxcat
	//2018�N12��8��
	
	//Please note that the following code is only for "Dynare for MATLAB".
	//You can use/rewrite this code without my permission.
	//I will not be responsible for any damage/liability induced by running this code.
	/////////////////////////////////////////////////////////////////////////////////////
	//����`RBC���f����Dynare�R�[�h
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
	model;
	c^(-theta)=beta*(1+r(+1)-delta)*c(+1)^(-theta);%����̃I�C���[������
	n^phi=w*c^(-theta);%�J���̍œK������
	y=z*(k(-1)^alpha)*(n^(1-alpha));%���Y�֐�
	r=z*alpha*(n/k(-1))^(1-alpha);%���{�̃����^����
	w=z*(1-alpha)*(k(-1)/n)^alpha;%������
	k=(1-delta)*k(-1)+i;%���{�̐��ڎ�
	y=c+i;%���s��̋ύt����
	log(z)=rho*log(z(-1))+e;%�Z�p�V���b�N
	end;
	
	//�p�[�g4�@����Ԓl�̌v�Z
	initval;%�����l�̐ݒ�
	c=1;
	n=1;
	y=3;
	r=0.05;
	w=1;
	k=10;
	i=1;
	z=1;
	end;
	steady;%����Ԃ̌v�Z
	check;%BK�����𖞂������m�F
	
	//�p�[�g5�@�C���p���X�����̌v�Z
	shocks; %�V���b�N�̐ݒ�
	var e; %�V���b�N��^����O���ϐ�
	periods 1;%�V���b�N��^������@//periods 1:10�Ȃǂ��\
	values 0.01;%�V���b�N�̑傫��
	end;
	
	simul(periods=100);%�V�~�����[�V�����̎��s�Ƃ��̊���
	
	//�p�[�g6�@����Ԃ���̘������̌v�Z�ƃv���b�g
	ctild=(c-c(1))/c(1);%����Ԃ���̘������̒�`
	ntild=(n-n(1))/n(1);
	ytild=(y-y(1))/y(1);
	rtild=(r-r(1))/r(1);
	wtild=(w-w(1))/w(1);
	ktild=(k-k(1))/k(1);
	itild=(i-i(1))/i(1);
	ztild=(z-z(1))/z(1);
	
	figure(1)
	plot(0:101,ytild,0:101,ctild,0:101,itild,0:101,ntild);
	%�����O���t��Ƀv���b�g�i0�`101���j
	
	figure(2)
	subplot(2,2,1)
	plot(0:100, ytild(1:101)); title('GDP')
	subplot(2,2,2)
	plot(0:100, ctild(1:101)); title('����')
	subplot(2,2,3)
	plot(0:100, itild(1:101)); title('����')
	subplot(2,2,4)
	plot(0:100, ntild(1:101)); title('�J��')
	
	figure(3)
	subplot(2,2,1)
	plot(0:100, wtild(1:101)); title('����')
	subplot(2,2,2)
	plot(0:100, ktild(1:101)); title('���{�~��')
	subplot(2,2,3)
	plot(0:100, rtild(1:101)); title('���{�̃����^���R�X�g')
	subplot(2,2,4)
	plot(0:100, ztild(1:101)); title('�Z�p����')
	
	//���ʂ�CSV�t�@�C���ɏ����o��
	csvwrite('rbc_rslt.csv',[ytild, ctild, itild, ntild, wtild, ktild, rtild, ztild]);