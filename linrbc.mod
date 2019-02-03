////////////////////////////////////////////////////////////////////////////////////
	//線形RBCモデル(カリブレーション)
	//Boxcat(東京大学経済学部４年) 
	//2018年12月8日
	
	//Please note that the following code is only for "Dynare for MATLAB".
	//You can use/rewrite this code without my permission.
	//I will not be responsible for any damage/liability induced by running this code.
	/////////////////////////////////////////////////////////////////////////////////////
//RBCモデルのDynareコード
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
model(linear);

%定常状態値で求めた数値
#rbar=(1+beta*delta-beta)/theta;%r/(1+r-delta);
#iybar=(delta*alpha*beta)/(1+beta*delta-beta);
#cybar=1-iybar;

y=z+alpha*k(-1)+(1-alpha)*n;%生産関数
r=z+(alpha-1)*k(-1)+(1-alpha)*n;%資本のレンタル料
w=z+alpha*k(-1)+alpha*n;%労働需要の最適化
c=c(+1)-rbar*r(+1);%消費のオイラー方程式
n=(1/phi)*w-(theta/phi)*c;%最適労働供給条件
y=cybar*c+iybar*i;%財市場の均衡条件
k=(1-delta)*k(-1)+delta*i;%資本の推移式
z=rho*z(-1)+e;%技術ショック
end;

//パート4　モデルのチェック
check;
resid(1);
steady;

//パート5　ショック
shocks;
var e = 0.01;
end;

//パート6　シミュレーション(100期)
stoch_simul(irf=100)y c i n w k r z;