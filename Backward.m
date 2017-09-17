function lmd = Backward( x, u, T, dv, q, sf, a, b, len )
    dt = T / dv;
    
    lmd = zeros( len.lmd * dv, 1 );
    lmd((1:len.lmd)+len.lmd*(dv-1)) = dPhidx( x((1:len.x)+len.x*(dv-1)), sf );
    
    for cnt = dv-1:-1:1
        lmd((1:len.lmd)+len.lmd*(cnt-1)) = lmd((1:len.lmd)+len.lmd*(cnt)) + dHdx( x((1:len.x)+len.x*(cnt)), lmd((1:len.lmd)+len.lmd*(cnt)), q, a, b ) * dt;
    end
end