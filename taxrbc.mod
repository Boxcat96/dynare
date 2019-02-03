////////////////////////////////////////////////////////////////////////////////////
//税制導入RBCモデル(non-linear)
//Boxcat
//2018年12月11日

//Please note that the following code is only for "Dynare for MATLAB".
//You can use/rewrite this code without my permission.
//I will not be responsible for any damage/liability induced by running this code.
/////////////////////////////////////////////////////////////////////////////////////

//変数の定義
var y c i n k rk w z g lambda; %内生変数
varexo tauc tauk taun ez;      %外生変数

//パラメータの定義
parameters beta theta varphi alpha delta rhoz;
parameters yss css iss nss kss rkss wss zss lambdass;

beta = 0.996;   %割引率
theta=1.5;      %消費の代替弾力性の逆数
varphi = 2.0;   %労働供給の代替弾力性の逆数
alpha = 0.33;   %資本分配率
delta = 0.04;   %資本減耗率
rhoz = 0.9;     %技術ショックの持続性

//モデル(非線形)
model;
c^(-theta) = (1+tauc)*lambda;       %消費のオイラー方程式
n^varphi = (1-taun)*w*lambda;       %労働の最適化条件
lambda = beta * lambda(+1)*(rk(+1)*(1-tauk)+1-delta); %ラグランジュ乗数
k = (1-delta)*k(-1) + i;            %資本の推移式
y = z*(k(-1)^alpha)*n^(1-alpha);    %生産関数
rk = z*(alpha)*(n/k(-1))^(1-alpha); %資本のレンタルコスト
w = z*(1-alpha)*(k(-1)/n)^(alpha);  %賃金率
y = c + i + g;                      %財市場の均衡条件
g = tauc*c + tauk*rk*k(-1) + taun*w*n;  %税収
log(z) = rhoz * log(z(-1)) + ez;    %技術ショック
end;

//定常状態値の定義(ss は steady state の意味です)
css=1;
nss=1;
yss=3;
rkss=0.05;
wss=1;
kss=10;
iss=1;
zss = 1;
lambdass = 1;

//初期値の設定
initval;
tauc = 0.05;  %初期時点での消費税
tauk = 0.2;   %初期時点での資本税
taun = 0.2;   %初期時点での所得税
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
steady;%初期の定常状態の計算

//終端値の設定
endval;
tauc = 0.1;   %終端時点での消費税
tauk = 0.2;   %終端時点での資本税
taun = 0.2;   %終端時点での所得税
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
steady;%ショック後の定常状態の計算

check;%BK条件を満たすか確認

shocks;
var tauc;         %ショックを与える外生変数
periods 0:5, 5:10;     %ショックの変更区間
values 0, 0.08;         %ショックの大きさ(0ならアナウンスメントだけ)
end;

simul(periods=100); %ショック(期間)

//定常状態からの乖離率の計算
ctild=(c-c(1))/c(1);
ntild=(n-n(1))/n(1);
ytild=(y-y(1))/y(1);
rktild=(rk-rk(1))/rk(1);
wtild=(w-w(1))/w(1);
ktild=(k-k(1))/k(1);
itild=(i-i(1))/i(1);
gtild=(g-g(1))/g(1);

//図の描画(定常状態からの乖離率のインパルス反応)
figure(1)
subplot(2,2,1)
plot(0:100, ytild(1:101)); title('GDP')
subplot(2,2,2)
plot(0:100, ctild(1:101)); title('消費')
subplot(2,2,3)
plot(0:100, itild(1:101)); title('投資')
subplot(2,2,4)
plot(0:100, wtild(1:101)); title('賃金')

figure(2)
subplot(2,2,1)
plot(0:100, ntild(1:101)); title('労働')
subplot(2,2,2)
plot(0:100, rktild(1:101)); title('資本のレンタルコスト')
subplot(2,2,3)
plot(0:100, ktild(1:101)); title('資本蓄積')
subplot(2,2,4)
plot(0:100, gtild(1:101)); title('税収')

//close all;