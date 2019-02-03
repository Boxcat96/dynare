////////////////////////////////////////////////////////////////////////////////////
	//非線形RBCモデル(カリブレーション)
	//Boxcat
	//2018年12月8日
	
	//Please note that the following code is only for "Dynare for MATLAB".
	//You can use/rewrite this code without my permission.
	//I will not be responsible for any damage/liability induced by running this code.
	/////////////////////////////////////////////////////////////////////////////////////
	//非線形RBCモデルのDynareコード
	//パート１　変数の定義
	var c n y r w k i z;%内生変数
	varexo e;%外生変数
	
	//パート２　パラメータの設定
	parameters alpha beta theta phi delta rho;%パラメータを定義
	//パラメーターに値を代入
	alpha=0.33;%資本分配率
	beta=0.996;%割引率
	theta=1.5;%消費の代替弾力性の逆数
	phi=2;%労働供給の代替弾力性の逆数
	delta=0.04;%資本減耗率
	rho=0.9;%技術進歩の持続性
	
	//パート3モデルの記述
	model;
	c^(-theta)=beta*(1+r(+1)-delta)*c(+1)^(-theta);%消費のオイラー方程式
	n^phi=w*c^(-theta);%労働の最適化条件
	y=z*(k(-1)^alpha)*(n^(1-alpha));%生産関数
	r=z*alpha*(n/k(-1))^(1-alpha);%資本のレンタル料
	w=z*(1-alpha)*(k(-1)/n)^alpha;%賃金率
	k=(1-delta)*k(-1)+i;%資本の推移式
	y=c+i;%財市場の均衡条件
	log(z)=rho*log(z(-1))+e;%技術ショック
	end;
	
	//パート4　定常状態値の計算
	initval;%初期値の設定
	c=1;
	n=1;
	y=3;
	r=0.05;
	w=1;
	k=10;
	i=1;
	z=1;
	end;
	steady;%定常状態の計算
	check;%BK条件を満たすか確認
	
	//パート5　インパルス反応の計算
	shocks; %ショックの設定
	var e; %ショックを与える外生変数
	periods 1;%ショックを与える期　//periods 1:10なども可能
	values 0.01;%ショックの大きさ
	end;
	
	simul(periods=100);%シミュレーションの実行とその期間
	
	//パート6　定常状態からの乖離率の計算とプロット
	ctild=(c-c(1))/c(1);%定常状態からの乖離率の定義
	ntild=(n-n(1))/n(1);
	ytild=(y-y(1))/y(1);
	rtild=(r-r(1))/r(1);
	wtild=(w-w(1))/w(1);
	ktild=(k-k(1))/k(1);
	itild=(i-i(1))/i(1);
	ztild=(z-z(1))/z(1);
	
	figure(1)
	plot(0:101,ytild,0:101,ctild,0:101,itild,0:101,ntild);
	%同じグラフ上にプロット（0～101期）
	
	figure(2)
	subplot(2,2,1)
	plot(0:100, ytild(1:101)); title('GDP')
	subplot(2,2,2)
	plot(0:100, ctild(1:101)); title('消費')
	subplot(2,2,3)
	plot(0:100, itild(1:101)); title('投資')
	subplot(2,2,4)
	plot(0:100, ntild(1:101)); title('労働')
	
	figure(3)
	subplot(2,2,1)
	plot(0:100, wtild(1:101)); title('賃金')
	subplot(2,2,2)
	plot(0:100, ktild(1:101)); title('資本蓄積')
	subplot(2,2,3)
	plot(0:100, rtild(1:101)); title('資本のレンタルコスト')
	subplot(2,2,4)
	plot(0:100, ztild(1:101)); title('技術水準')
	
	//結果をCSVファイルに書き出し
	csvwrite('rbc_rslt.csv',[ytild, ctild, itild, ntild, wtild, ktild, rtild, ztild]);