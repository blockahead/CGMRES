%% 函数Hの状態微分
% x        : [ x;dx ]       （位置，速度）
% lmd      : [ lmd1;lmd2 ]  （随伴変数1，随伴変数2）
% q        : [ q1;q2 ]      （位置の重み，速度の重み）
% a        : a              （システム係数）
% b        : b              （システム係数）

function Hx = dHdx( x, lmd, sys, cgmres )
    Hx = [ ...
        sys.a * lmd(2) + cgmres.q(1) * x(1); ...
        sys.b * lmd(2) + cgmres.q(2) * x(2) + lmd(1); ...
        ];
end