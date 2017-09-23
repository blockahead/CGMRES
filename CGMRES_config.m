close all;
clear;


% �V�X�e����`
sys.a = -1;     % �V�X�e���ϐ�
sys.b = -1;     % �V�X�e���ϐ�

% �V�~�����[�V������`
tsim = 20; % �V�~�����[�V�������� (s)


% C/GMRES�̃R���g���[����`
dSamplingPeriod = 0.001;                % MATLAB Function�̃T���v�����O���� (s)
cgmres.ht = 0.001;                      % �O�i�����ߎ��̎��ԕ� (s)
cgmres.zeta = 1000.0;                   % ����ʂ̈��艻�Q�C�� (-)
cgmres.tf = 1.0;                        % �\�����Ԃ̍ŏI�l (s)
cgmres.alpha = 0.5;                     % �\�����Ԃ̏㏸���x�Q�C�� (-)
cgmres.dv = 5;                          % �\�����Ԃ̕����� (-) �i�]�������ɂ���ĕ]������|�C���g�̐��j

cgmres.x0 = [2;0];                      % �R���g���[���ɗ^���鏉�����
cgmres.u0 = [0.01;0.9;0.03];            % �R���g���[���ɗ^���鏉�������

cgmres.len_x = length( cgmres.x0 );     % ��Ԃ̐�
cgmres.len_u = length( cgmres.u0 );     % ����ʂ̐�
cgmres.len_lmd = cgmres.len_x;          % �����ϐ��̐�

cgmres.q = [ 1;10 ];                    % ��Ԃɑ΂���d��
cgmres.r = [ 1;0.01 ];                  % ����ʂɑ΂���d��
cgmres.sf = [ 1;10 ];                   % �\�����Ԃ̍ŏI��Ԃɑ΂���d��

cgmres.umax = 1;                        % ���͏���i�����̓[���ɐݒ肵�Ă���j

% C/GMRES�̃R���g���[���p�v�Z
buff = c2d( ss( tf( [ cgmres.tf ], [ (1/cgmres.alpha), 1 ] ) ), dSamplingPeriod );
cgmres.T_outGain = buff.a;               % �\�����ԍ���������
cgmres.T_inGain = buff.b;              % �\�����ԍ���������
clearvars buff;

% �V�~�����[�V�����̎��s�Ǝ��Ԍv��
tic;
sim( 'CGMRES' );
toc