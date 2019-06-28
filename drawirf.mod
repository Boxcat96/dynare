s////////////////////////////////////////////////////////////////////////////////////
	//線形RBCモデル(カリブレーション)
	//Boxcat 
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

 
set_param_value('rho',0.5)
stoch_simul(irf=7,order=1,irf_plot_threshold=0,nograph,noprint) y r c i ; 
irf_0_trend=oo_.irfs;
set_param_value('rho',0.7)
stoch_simul(irf=7,order=1,irf_plot_threshold=0,nograph,noprint) y r c i ; 
irf_2_trend=oo_.irfs;
set_param_value('rho',0.8)
stoch_simul(irf=7,order=1,irf_plot_threshold=0,nograph,noprint) y r c i ; 
irf_4_trend=oo_.irfs;
set_param_value('rho',0.9)
stoch_simul(irf=7,order=1,irf_plot_threshold=0,nograph,noprint) y r c i ; 
irf_6_trend=oo_.irfs;

figure('Name',' IRF ')
plot(0:options_.irf,[0 irf_0_trend.y_e],'k-',0:options_.irf,[0 irf_2_trend.y_e],'b--',0:options_.irf,[0 irf_4_trend.y_e],'r-.',0:options_.irf,[0 irf_6_trend.y_e],'*-')
