%% ����ʂ̕ω��ʂ̌v�Z�iC/GMRES�@�j
function du_new = GMRES( x_current, du, u, T, sys, cgmres )
    
    k = cgmres.len_u * cgmres.dv; % GMRES�̍ő�J��Ԃ��񐔁i���U����2�_���E�l���̗v�f���Ɠ������j

    
    % GMRES�@�ŉ����ׂ����̐ݒ�
    
    dx = Func( x_current, u, sys ); % ��ԕω��ʂ̑O�i�����ߎ�
    
    Fxt = CalcF( x_current + ( dx * cgmres.ht ), u, T, sys, cgmres ); % ����F�̏�Ԕ���
    
    F = CalcF( x_current, u, T, sys, cgmres ); % ����F
    
    Right = - cgmres.zeta * F - ( ( Fxt - F ) / cgmres.ht ); % ���m�����E�ӂɂ܂Ƃ߂�
    
    Fuxt = CalcF( x_current + ( dx * cgmres.ht ), u + ( du * cgmres.ht ), T, sys, cgmres ); % ����F�̓��́C��Ԕ���
    
    Left = ( ( Fuxt - Fxt ) / cgmres.ht ); % ���m�������ӂɂ܂Ƃ߂�
    
    
    % GMRES�@�ł̌v�Z
    
    r0 = Right - Left; % �����c��
    
    du_new = zeros( cgmres.len_u * cgmres.dv, 1 );
    
    v = zeros( cgmres.len_u * cgmres.dv, k + 1 ); 
    v(:,1) = r0 / norm( r0 );

    h = zeros( k + 1 );

    y = zeros( cgmres.len_u * cgmres.dv, 1 );
    y_pre = zeros( cgmres.len_u * cgmres.dv, 1 );
    
    e = zeros( k + 1, 1 );
    e(1) = 1;
    
    for cnt = 1:k
        Fuxt = CalcF( x_current + ( dx * cgmres.ht ), u + ( v(:,cnt) * cgmres.ht ), T, sys, cgmres );
        Av = ( ( Fuxt - Fxt ) / cgmres.ht );

        Sum = zeros( cgmres.len_u * cgmres.dv, 1 );
    
        for cnt2 = 1:cnt
            h(cnt2,cnt) = Av' * ( v(:,cnt2) );
            Sum = Sum + h(cnt2,cnt) * v(:,cnt2);
        end

        v_est = Av - Sum;

        h(cnt+1,cnt) = norm( v_est );

        v(:,cnt+1) = v_est / h(cnt+1,cnt);

        y(1:cnt) = h(1:(cnt+1),1:cnt) \ ( norm( r0 ) * e(1:(cnt+1)) );
        
        % �덷���������Ȃ�����ł��؂�
        if ( norm( norm( r0 ) * e(1:(cnt+1)) - h(1:(cnt+1),1:cnt) * y(1:cnt) ) < 0.001 )
            du_new = du + v(:,1:(cnt-1)) * y_pre(1:(cnt-1));
            break;
        end
        
        y_pre = y;
        
    end
end