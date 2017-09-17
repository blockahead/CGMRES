% x        : [ x;dx ]       （位置，速度）
% u        : [ u;v;mu ]     （入力，疑似変数，ラグランジュ乗数）
% lmd      : [ lmd1;lmd2 ]  （随伴変数1，随伴変数2）
% r        : [ r1;r2 ]      （位置の重み，速度の重み）
% umax     : umax           （入力の最大値）
% a        : a              （システム係数）
% b        : b              （システム係数）

function Hu = dHdu( x, u, lmd, r, umax, a, b )
    Hu = [ ...
        r(1) * u(1) + b * lmd(2) * x(2) + u(3) * ( 2 * u(1) - umax );
        2 * u(2) * u(3) - r(2);
        ( u(1) - ( umax / 2 ) )^2 + u(2)^2 - ( umax / 2 )^2;
        ];
end