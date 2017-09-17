function F = CalcF( x_current, u, T, dv, q, r, sf, a, b, umax, len )
    x = Forward( x_current, u, T, dv, a, b, len );
    lmd = Backward( x, u, T, dv, q, sf, a, b, len );

    F = zeros( len.u * dv, 1 );

    for cnt = 1:dv
        F((1:len.u)+len.u*(cnt-1)) = dHdu(x((1:len.x)+len.x*(cnt-1)), u((1:len.u)+len.u*(cnt-1)), lmd((1:len.lmd)+len.lmd*(cnt-1)), r, umax, a, b );
    end
end