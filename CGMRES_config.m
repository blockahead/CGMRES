close all;
clear;


%% システム定義
sys.a = -1;     % システム変数
sys.b = -1;     % システム変数

%% シミュレーション定義
tsim = 20; % シミュレーション時間 (s)


%% C/GMRESのコントローラ定義
dSamplingPeriod = 0.001;                % MATLAB Functionのサンプリング周期 (s)
cgmres.ht = 0.001;                      % 前進差分近似の時間幅 (s)
cgmres.zeta = 1000.0;                   % 操作量の安定化ゲイン (-)
cgmres.tf = 1.0;                        % 予測時間の最終値 (s)
cgmres.alpha = 0.5;                     % 予測時間の上昇速度ゲイン (-)
cgmres.dv = 5;                          % 予測時間の分割数 (-) （評価函数によって評価するポイントの数）

cgmres.x0 = [2;0];                      % コントローラに与える初期状態
cgmres.u0 = [0.01;0.9;0.03];            % コントローラに与える初期操作量



cgmres.q = [ 1;10 ];                    % 状態に対する重み
cgmres.r = [ 1;0.01 ];                  % 操作量に対する重み
cgmres.sf = [ 1;10 ];                   % 予測時間の最終状態に対する重み

cgmres.umax = 1;                        % 入力上限（下限はゼロに設定している）

%% C/GMRESのコントローラ用計算
% 予測時間Tの計算用定数
buff = c2d( ss( tf( [ cgmres.tf ], [ (1/cgmres.alpha), 1 ] ) ), dSamplingPeriod );
cgmres.T_outGain = buff.a;              % 予測時間差分方程式
cgmres.T_inGain = buff.b;               % 予測時間差分方程式
clearvars buff;

% 初期入力値の計算（Newton法）
lmd0 = dPhidx( cgmres.x0, cgmres );
u0 = [1;2;3;]; % Newton法の初期値

for cnt = 1:20
    cgmres.u0 = cgmres.u0 - ddHddu( cgmres.x0, cgmres.u0, lmd0, sys, cgmres ) \ dHdu( cgmres.x0, cgmres.u0, lmd0, sys, cgmres );
end

cgmres.len_x = length( cgmres.x0 );     % 状態の数
cgmres.len_u = length( cgmres.u0 );     % 操作量の数
cgmres.len_lmd = cgmres.len_x;          % 随伴変数の数

%% シミュレーションの実行と時間計測
tic;
sim( 'CGMRES' );
toc