%% ����H�̓��͔���
% x        : [ x;dx ]       �i�ʒu�C���x�j
% u        : [ u;v;mu ]     �i����ʁC�_�~�[����ʁC���O�����W���搔�j
% lmd      : [ lmd1;lmd2 ]  �i�����ϐ�1�C�����ϐ�2�j
% r        : [ r1;r2 ]      �i����ʂ̏d�݁C�_�~�[����ʂ̏d�݁j
% umax     : umax           �i����ʂ̍ő�l�j
% a        : a              �i�V�X�e���W���j
% b        : b              �i�V�X�e���W���j

function Hu = dHdu( x, u, lmd, sys, cgmres )
    Hu = [ ...
        cgmres.r(1) * u(1) + sys.b * lmd(2) * x(2) + u(3) * ( 2 * u(1) - cgmres.umax );
        2 * u(2) * u(3) - cgmres.r(2);
        ( u(1) - ( cgmres.umax / 2 ) )^2 + u(2)^2 - ( cgmres.umax / 2 )^2;
        ];
end