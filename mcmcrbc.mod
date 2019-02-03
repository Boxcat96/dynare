////////////////////////////////////////////////////////////////////////////////////
//RBCモデル(ベイズ推計) Boxcat
//2018年12月8日

//Please note that the following code is only for "Dynare for MATLAB".
//You can use/rewrite this code without my permission.
//I will not be responsible for any damage/liability induced by running this code.
/////////////////////////////////////////////////////////////////////////////////////
//パート１　変数の定義
var c n y r w k i z y_obs c_obs i_obs w_obs n_obs;%内生変数
varexo e uy uc ui uw un;%外生変数

//パート２　パラメータの設定
parameters alpha beta theta phi delta rho;%パラメータを定義
//パラメーターに値を代入
beta=0.996;%割引率
delta=0.04;%資本減耗率

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

y_obs=y+uy; %GDPの観測方程式
c_obs=c+uc; %消費の観測方程式
i_obs=i+ui; %投資の観測方程式
w_obs=w+uw; %賃金率の観測方程式
n_obs=n+un; %労働量の観測方程式
end;

//パート4　推測するパラメータの事前分布
estimated_params;
alpha, beta_pdf, 0.33, 0.1;%資本分配率
theta, gamma_pdf, 1.5, 0.1;%消費の代替弾力性の逆数
phi, gamma_pdf, 2, 0.1;%労働供給の代替弾力性の逆数
rho, beta_pdf, 0.85, 0.1;%技術ショックの持続性
stderr e, inv_gamma_pdf, 0.1, inf;
stderr uy, inv_gamma_pdf, 0.1, inf;
stderr uc, inv_gamma_pdf, 0.1, inf;
stderr ui, inv_gamma_pdf, 0.1, inf;
stderr uw, inv_gamma_pdf, 0.1, inf;
stderr un, inv_gamma_pdf, 0.1, inf;
end;
%チェックコマンドは使ってはならない。

%パート5　シミュレーション
varobs y_obs c_obs i_obs w_obs n_obs;%読み取る観測データ

estimation(datafile=rbcdata, mode_check, mh_replic=100000, mh_nblocks=2,
mh_drop=0.5, mh_jscale=0.73, bayesian_irf);

stoch_simul;