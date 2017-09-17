%% ����ʂ̕ω��ʂ̌v�Z�iC/GMRES�@�j
function du_new = GMRES( x_current, du, u, T, dv, q, r, sf, zeta, a, b, umax, ht, len )
    
    k = len.u * dv; % GMRES�̍ő�J��Ԃ��񐔁i���U����2�_���E�l���̗v�f���Ɠ������j

    
    % GMRES�@�ŉ����ׂ����̐ݒ�
    
    dx = Func( x_current, u, a, b ); % ��ԕω��ʂ̑O�i�����ߎ�
    
    Fxt = CalcF( x_current + ( dx * ht ), u, T, dv, q, r, sf, a, b, umax, len ); % ����F�̏�Ԕ���
    
    F = CalcF( x_current, u, T, dv, q, r, sf, a, b, umax, len ); % ����F
    
    Right = - zeta * F - ( ( Fxt - F ) / ht ); % ���m�����E�ӂɂ܂Ƃ߂�
    
    Fuxt = CalcF( x_current + ( dx * ht ), u + ( du * ht ), T, dv, q, r, sf, a, b, umax, len ); % ����F�̓��́C��Ԕ���
    
    Left = ( ( Fuxt - Fxt ) / ht ); % ���m�������ӂɂ܂Ƃ߂�
    
    
    % GMRES�@�ł̌v�Z
    
    r0 = Right - Left; % �����c��
    
    du_new = zeros( len.u * dv, 1 );
    
    v = zeros( len.u * dv, k + 1 ); 
    v(:,1) = r0 / norm( r0 );

    h = zeros( k + 1 );

    y = zeros( len.u * dv, 1 );
    y_pre = zeros( len.u * dv, 1 );
    
    e = zeros( k + 1, 1 );
    e(1) = 1;
    
    for cnt = 1:k
        Fuxt = CalcF( x_current + ( dx * ht ), u + ( v(:,cnt) * ht ), T, dv, q, r, sf, a, b, umax, len );
        Av = ( ( Fuxt - Fxt ) / ht );

        Sum = zeros( len.u * dv, 1 );
    
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