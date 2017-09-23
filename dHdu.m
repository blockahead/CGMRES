%% 函数Hの入力微分
% x        : [ x;dx ]       （位置，速度）
% u        : [ u;v;mu ]     （操作量，ダミー操作量，ラグランジュ乗数）
% lmd      : [ lmd1;lmd2 ]  （随伴変数1，随伴変数2）
% r        : [ r1;r2 ]      （操作量の重み，ダミー操作量の重み）
% umax     : umax           （操作量の最大値）
% a        : a              （システム係数）
% b        : b              （システム係数）

function Hu = dHdu( x, u, lmd, sys, cgmres )
    Hu = [ ...
        cgmres.r(1) * u(1) + sys.b * lmd(2) * x(2) + u(3) * ( 2 * u(1) - cgmres.umax );
        2 * u(2) * u(3) - cgmres.r(2);
        ( u(1) - ( cgmres.umax / 2 ) )^2 + u(2)^2 - ( cgmres.umax / 2 )^2;
        ];
end