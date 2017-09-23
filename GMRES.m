%% 操作量の変化量の計算（C/GMRES法）
function du_new = GMRES( x_current, du, u, T, sys, cgmres )
    
    k = cgmres.len_u * cgmres.dv; % GMRESの最大繰り返し回数（離散時間2点境界値問題の要素数と等しい）

    
    % GMRES法で解くべき問題の設定
    
    dx = Func( x_current, u, sys ); % 状態変化量の前進差分近似
    
    Fxt = CalcF( x_current + ( dx * cgmres.ht ), u, T, sys, cgmres ); % 函数Fの状態微分
    
    F = CalcF( x_current, u, T, sys, cgmres ); % 函数F
    
    Right = - cgmres.zeta * F - ( ( Fxt - F ) / cgmres.ht ); % 既知項を右辺にまとめる
    
    Fuxt = CalcF( x_current + ( dx * cgmres.ht ), u + ( du * cgmres.ht ), T, sys, cgmres ); % 函数Fの入力，状態微分
    
    Left = ( ( Fuxt - Fxt ) / cgmres.ht ); % 未知項を左辺にまとめる
    
    
    % GMRES法での計算
    
    r0 = Right - Left; % 初期残差
    
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
        
        % 誤差が小さくなったら打ち切り
        if ( norm( norm( r0 ) * e(1:(cnt+1)) - h(1:(cnt+1),1:cnt) * y(1:cnt) ) < 0.001 )
            du_new = du + v(:,1:(cnt-1)) * y_pre(1:(cnt-1));
            break;
        end
        
        y_pre = y;
        
    end
end