%% ����H�̏�Ԕ���
% x        : [ x;dx ]       �i�ʒu�C���x�j
% lmd      : [ lmd1;lmd2 ]  �i�����ϐ�1�C�����ϐ�2�j
% sys      : a              �i�V�X�e���W���j
% sys      : b              �i�V�X�e���W���j
% cgmres   : [ q1;q2 ]      �i�ʒu�̏d�݁C���x�̏d�݁j

function Hx = dHdx( x, lmd, u, sys, cgmres )
    Hx = [ ...
        sys.a * lmd(2) + cgmres.q(1) * x(1); ...
        sys.b * lmd(2) * u(1) + cgmres.q(2) * x(2) + lmd(1); ...
        ];
end