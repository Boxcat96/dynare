////////////////////////////////////////////////////////////////////////////////////
//変更版CEE(Christiano, Eichenbaum and Evans, 2005)モデル
//中規模マクロDSGEモデル(カリブレーション)
//2018年12月7日

//Please note that the following code is only for "Dynare for MATLAB".
//You can use/rewrite this code without my permission.
//I will not be responsible for any damage/liability induced by running this code.
/////////////////////////////////////////////////////////////////////////////////////

//パート１：変数の定義///////////////////////////////////

//内生変数
var y c co cr lambda i q k kg n w rk b r pi g gi tau z zeta ucc uii unn;

//外生変数
varexo eg ez et em eq ec ei egi ew en;

//パラメータ
parameters beta theta varphi h omega alpha nu deltap deltag lw rw zw;
parameters kappa eta psi phi phih phir phipi phiy phit phiyt phibt phiig;
parameters rhog rhoz rhoc rhoi rhogi rhon bybar gybar;

//パート２：カリブレートするパラメータの定義(江口2011)//////////////////////////////////

h    = 0.7742 ; %消費者の習慣形成の度合い
beta = 0.996; %割引率
theta = 1.5; %消費の代替弾力性
varphi = 2; %労働供給の代替弾力性
alpha = 0.33; %資本分配率
omega   =     0.4230   ;   %非リカード家計の割合
deltap = 0.06; %民間資本減耗率
deltag = 0.04; %社会資本減耗率
eta = 0.75; %Calvoパラメーター
kappa = 7; %投資の調整コスト
psi = 11; %マークアップ率
zw  = 0.3760      ;  %賃金を最適化できない確率
lw  = 0.2; %賃金マークアップショック
rw  = 0.2; %賃金の過去のインフレ率に対する連動度合い
nu      =     0.2520    ;   %社会資本の生産力効果
phi      =     0.6509    ;   %金融政策ショックの持続性
phih    = 0.1231  ; %物価の過去のインフレ率に対する連動度合い(ハイブリッドNKPC)
phipi     =    1.7423     ;  %テイラールールのインフレ率の項
phir      =   0.0835       ; %テイラールールの名目金利の持続性
phiy      =   0.3272       ; %テイラールールのGDP項
phit      =   0.6382       ; %課税ルールの前期項
phiig      = 0.0050        ;  %公共投資ショックの前期GDP項
phibt     =   0.1794       ; %課税ルールの国債の項
phiyt      =  0.0070      ; %課税ルールのGDPの項
rhog   =      0.5523      ; %公共消費ショックの持続性
rhogi        =0.4166      ; %公共投資ショックの持続性
rhoz    =     0.8765       ; %技術ショックの持続性
rhoc         =0.4737      ; %選好ショックの持続性
rhoi         =0.49      ; %投資ショックの持続性
rhon         =0.51       ; %労働供給ショックの持続性

bybar = 2.35; %定常における国債の対GDP比
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
q=pi(+1)-r+((1-deltap)/(1-deltap+rkbar))*q(+1)+(rkbar/(1-deltap
+rkbar))*rk(+1)+eq; %トービンのq
//pi=beta*pi(+1)+(((1-eta)*(1-beta*eta))/(eta))*((1-alpha)*w+alpha*rk-
//nu*kg(-1)-z); %標準NKPC
//ハイブリッドNKPC///////////////////////////////////////
pi=(beta/(1+beta*phih))*pi(+1)+(phih/
1+beta*phih)*pi(-1)+(1/1+beta*phih)*(((1-eta)*(1-beta*eta))/(eta))*((1-
alpha)*w+alpha*rk-nu*kg(-1)-z); %ハイブリッドNKPC
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
unn=rhon*unn(-1)+en; %労働供給ショック
end;

//パート４：ショック、シミュレーション/////////////////////////////////////
shocks;
var et = 1; %減税ショック
var em = 1; %金融緩和ショック
var ez = 1; %技術ショック
var eg = 1; %政府消費ショック
var eq = 1; %証券プレミアムショック
var ec = 1; %選好ショック
var ei = 1; %投資ショック
var ew = 1; %賃金マークアップショック
var en = 1; %労働供給ショック
var egi = 1; %公共投資ショック
end;

check; resid(1); steady; %BK条件の確認、定常状態値 

stoch_simul(irf=40) y c i n k b w r pi ; %確率的シミュレーション(40期)
//close all; %//を消すと自動的にグラフが消えるようになる(計算だけされる)。