%% ����H�̏�Ԕ���
% x        : [ x;dx ]       �i�ʒu�C���x�j
% lmd      : [ lmd1;lmd2 ]  �i�����ϐ�1�C�����ϐ�2�j
% q        : [ q1;q2 ]      �i�ʒu�̏d�݁C���x�̏d�݁j
% a        : a              �i�V�X�e���W���j
% b        : b              �i�V�X�e���W���j

function Hx = dHdx( x, lmd, sys, cgmres )
    Hx = [ ...
        sys.a * lmd(2) + cgmres.q(1) * x(1); ...
        sys.b * lmd(2) + cgmres.q(2) * x(2) + lmd(1); ...
        ];
end