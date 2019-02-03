//////////////////////////////////////////////////////////////////////////////
//変更版CEE(Christiano, Eichenbaum and Evans, 2005)モデル
//中規模マクロ計量DSGEモデル
//2018年12月7日

//Please note that the following code is only for "Dynare for MATLAB".
//You can use/rewrite this code without my permission.
//I will not be responsible for any damage/liability induced by running this code.
////////////////////////////////////////////////////////////////////////////

//パート１：変数の定義///////////////////////////////////

//内生変数
var y c co cr lambda i q n k kg w rk b r pi g gi tau z zeta ucc uii unn;
var y_obs c_obs i_obs g_obs pi_obs r_obs w_obs n_obs b_obs;

//外生変数
varexo ec ei ew en eg egi ez et em eq uy uc ui ug upi ur uw un ub;

//パラメータ
parameters beta theta varphi h omega alpha nu deltap deltag lw rw zw;
parameters kappa eta psi phi phih phir phipi phiy phit phiyt phibt phiig;
parameters rhog rhoz rhoc rhoi rhogi rhon bybar gybar;

//パート２：カリブレートするパラメータの定義(江口2011)//////////////////////////////////

beta = 0.996; %割引率
theta = 1.5; %消費の代替弾力性
varphi = 2; %労働供給の代替弾力性
alpha = 0.33; %資本分配率
deltap = 0.06; %民間資本減耗率
deltag = 0.04; %社会資本減耗率
lw=0.2; %賃金マークアップショック
rw=0.2; %賃金の過去のインフレ率に対する連動度合い
eta = 0.75; %Calvoパラメーター
kappa = 7; %投資の調整コスト
psi = 11; %マークアップ率

bybar = 4; %定常における国債の対GDP比
gybar = 0.08; %定常における公共投資の対GDP比

//パート３：モデル///////////////////////////////////////

model(linear);

//定常状態値の定義

# rbar = 1/beta; %国債粗利子率の定常状態値
# rkbar=rbar+deltap-1; %資本のレンタル料
# wbar=(((psi-1)/psi)*((((1-alpha)^(1-alpha))*(alpha^(alpha)))/
(rkbar^(alpha))))^(1/(1-alpha)); %定常における賃金率
# nybar=((((1-alpha)*rkbar)/(alpha*wbar)))^(alpha); %定常における労働量の対GDP比
# kybar=((((1-alpha)*rkbar)/(alpha*wbar)))^(alpha-1); %定常における資本の対GDP比
# iybar=deltap*kybar; %定常における消費の対GDP比
# cybar=1-iybar-gybar; %定常における投資の対GDP比

//線型モデル

(1-h)*(1-beta*h)*lambda = -theta*(co-h*co(-1))+(1-h)*ucc
+beta*h*(theta*(co(+1)-h*co)-(1-h)*ucc(+1)); %リカード的家計の消費
lambda = lambda(+1)+r-pi(+1);
cr=((wbar*nybar)/cybar)*(w+n)-(1/cybar)*tau; %非リカード的家計の消費
c=omega*cr+(1-omega)*co; %経済全体の消費
//実質賃金////////////////////////////////////////////////////
w=(beta/(1+beta))*w(+1)+(1/(1+beta))*w(-1)+(beta/(1+beta))*pi(+1)-
((1+beta*rw)/1+beta)*pi+
(rw/(1+beta))*pi(-1)-(1/(1+beta))*((lw*(1-beta*zw)*(1-zw))/(lw
+(1+lw)*varphi)*zw)*(w-varphi*n+lambda-unn)-ew;
/////////////////////////////////////////////////////////////
i=(1/(1+beta))*i(-1)+(beta/(1+beta))*i(+1)+(kappa/(1+beta))*q+(beta/
(1+beta))*(uii(+1)-uii); %投資関数
q=pi(+1)-r+((1-deltap)/(1-deltap+rkbar))*q(+1)+(rkbar/(1-deltap+rkbar))*rk(+1)
+eq; %トービンのq
//pi=beta*pi(+1)+(((1-eta)*(1-beta*eta))/(eta))*((1-alpha)*w+alpha*rk-
//nu*kg(-1)-z); %標準NKPC
//ハイブリッドNKPC///////////////////////////////////////
pi=(beta/(1+beta*phih))*pi(+1)+(phih/1+beta*phih)*pi(-1)+(1/1+beta*phih)*(((1-
eta)*(1-beta*eta))/(eta))*((1-alpha)*w+alpha*rk-nu*kg(-1)-z); %ハイブリッドNKPC
////////////////////////////////////////////////////////
k=(1-deltap)*k(-1)+deltap*i; %資本の推移式
y=z+alpha*k(-1)+(1-alpha)*n+nu*kg(-1); %生産関数
y=cybar*c+iybar*i+g +gi; %財市場の均衡条件
r=phir*r(-1)+(1-phir)*(phipi*((pi+pi(-1)+pi(-2)+pi(-3))/4)+phiy*y)+zeta; 
%テイラールール(金利スムージング)
n-k(-1)=rk-w; %費用最小化条件
b=rbar*b(-1)+bybar*rbar*r(-1)-rbar*bybar*pi+gi+g-tau; %政府の予算制約
kg=(1-deltag)*kg(-1)+deltag*(1/gybar)*g; %社会資本の蓄積式
tau=phit*tau(-1)+(1-phit)*(phiyt*y+phibt*b(-1))-et; %持続性のある税制ルール
gi=rhogi*gi(-1)+phiig*y(-1)+egi; %公共投資ショック
g=rhog*g(-1)+eg; %政府消費ショック
z=rhoz*z(-1)+ez; %技術ショック
zeta=phi*zeta(-1)-em; %金融政策ショック
ucc=rhoc*ucc(-1)+ec; %選好ショック
uii=rhoi*uii(-1)+ei; %投資ショック
unn=rhon*unn(-1)+en; %労働ショック

//観測方程式

y_obs=y+uy; %GDPの観測方程式
c_obs=c+uc; %消費の観測方程式
i_obs=i+ui; %投資の観測方程式
g_obs=g+gi+ug; %政府支出の観測方程式
pi_obs=pi+upi; %インフレ率の観測方程式
r_obs=r+ur; %名目利子率の観測方程式
w_obs=w+uw; %賃金率の観測方程式
n_obs=n+un; %労働量の観測方程式
b_obs=b+ub; %国債の観測方程式
end;

//パート４：推計するパラメータの事前分布(江口2011、廣瀬2012、鎌田2009など)///////////

estimated_params;

h, normal_pdf, 0.8, 0.15; %消費者の習慣形成の度合い
omega, beta_pdf, 0.3, 0.1; %非リカード家計の割合
nu, normal_pdf, 0.25, 0.1; %社会資本の生産力効果
zw, beta_pdf, 0.375, 0.1; %賃金を最適化できない確率
phi, normal_pdf, 0.22, 0.1; %金融政策ショックの持続性
phipi, normal_pdf, 1.775, 0.2; %テイラールールのインフレ率の項
phir, normal_pdf, 0.02, 0.1; %テイラールールの名目金利の持続性
phiy, normal_pdf, 0.1, 0.1; %テイラールールのGDP項
phih, beta_pdf, 0.25, 0.1; 
%物価の過去のインフレ率に対する連動度合い(ハイブリッドNKPC)
phit, beta_pdf, 0.8, 0.1; %課税ルールの前期項
phiig, normal_pdf, 0.1, 0.1; %公共投資ショックの前期GDP項
phibt, normal_pdf, 0.1, 0.1; %課税ルールの国債の項
phiyt, normal_pdf, 0.12, 0.1; %課税ルールのGDPの項
rhog, beta_pdf, 0.6, 0.1; %公共消費ショックの持続性
rhoz, beta_pdf, 0.9, 0.1; %技術ショックの持続性
rhoc, beta_pdf, 0.5, 0.2; %選好ショックの持続性
rhoi, beta_pdf, 0.5, 0.2; %投資ショックの持続性
rhogi, beta_pdf, 0.5, 0.2; %公共投資ショックの持続性
rhon, beta_pdf, 0.5, 0.2; %労働ショックの持続性
stderr et, inv_gamma_pdf, 0.1, inf; %減税ショック
stderr eg, inv_gamma_pdf, 0.1, inf; %政府消費ショック
stderr ez, inv_gamma_pdf, 0.1, inf; %技術ショック
stderr em, inv_gamma_pdf, 0.1, inf; %金融緩和ショック
stderr uy, inv_gamma_pdf, 0.1, inf; %以下、各内生変数の観測誤差
stderr uc, inv_gamma_pdf, 0.1, inf;
stderr ui, inv_gamma_pdf, 0.1, inf;
stderr ug, inv_gamma_pdf, 0.1, inf;
stderr upi, inv_gamma_pdf, 0.1, inf;
stderr ur, inv_gamma_pdf, 0.1, inf;
stderr uw, inv_gamma_pdf, 0.1, inf;
stderr un, inv_gamma_pdf, 0.1, inf;
stderr ub, inv_gamma_pdf, 0.1, inf;
end;

//パート５：推計///////////////////////////////////////////////////////////

//観測する変数の指定
varobs y_obs c_obs i_obs g_obs pi_obs r_obs w_obs n_obs b_obs;

//推計コマンド(カルマンスムージングあり)

//estimation(
//datafile=ch8data, %データファイル名(1983～1998、日本)
//mode_check, %モードチェック
//mh_replic=100000, %試行回数
//mh_nblocks=2, %チェーンの数
//mh_drop=0.5, mh_jscale=0.37, %受容率に影響
//bayesian_irf, %ベイズ推定 
//smoother %カルマンスムージング
//);

//推計コマンド(カルマンスムージングなし)//////////////////////////////
estimation(datafile=ch8data, mode_check, mh_replic=10000, mh_nblocks=2, 
mh_drop=0.5, mh_jscale=0.37, bayesian_irf);
////////////////////////////////////////////////////////////////////////


//合理的期待均衡解の導出
stoch_simul;

//ヒストリカル分解を行うコマンド
shock_decomposition(parameter_set=posterior_mean) y_obs;
