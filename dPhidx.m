%% �w�b�Z�s��̏�Ԕ���
% x        : [ x;dx ]       �i�ʒu�C���x�j
% sf       : [ sf1;sf2 ]    �i����T�ł̈ʒu�̏d�݁C����T�ł̑��x�̏d�݁j

function Phix = dPhidx( x, cgmres )
    Phix = [ x(1) * cgmres.sf(1);x(2) * cgmres.sf(2) ];
end