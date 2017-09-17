function du_new = GMRES( x_current, du, u, T, dv, q, r, sf, zeta, a, b, umax, ht, len )
    
    % GMRESの最大繰り返し回数（離散時間2点境界値問題の要素数と等しい）
    k = len.u * dv;

    % 前進差分近似
    dx = Func( x_current, u, a, b );
    
    Fxt = CalcF( x_current + ( dx * ht ), u, T, dv, q, r, sf, a, b, umax, len );
    F = CalcF( x_current, u, T, dv, q, r, sf, a, b, umax, len );
    Right = - zeta * F - ( ( Fxt - F ) / ht );
    
    Fuxt = CalcF( x_current + ( dx * ht ), u + ( du * ht ), T, dv, q, r, sf, a, b, umax, len );
    Left = ( ( Fuxt - Fxt ) / ht );
    
    r0 = Right - Left;
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
        
        % 打ち切り
        if ( norm( norm( r0 ) * e(1:(cnt+1)) - h(1:(cnt+1),1:cnt) * y(1:cnt) ) < 0.001 )
            du_new = du + v(:,1:(cnt-1)) * y_pre(1:(cnt-1));
            break;
        end
        
        y_pre = y;
        
    end
end