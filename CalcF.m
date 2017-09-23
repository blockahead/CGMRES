%% C/GMRES‚Å‚Ì”Ÿ”F‚ÌŒvZ
function F = CalcF( x_current, u, T, sys, cgmres )
    x = Forward( x_current, u, T, sys, cgmres );
    lmd = Backward( x, u, T, sys, cgmres );

    F = zeros( cgmres.len_u * cgmres.dv, 1 );

    for cnt = 1:cgmres.dv
        F((1:cgmres.len_u)+cgmres.len_u*(cnt-1)) = dHdu(x((1:cgmres.len_x)+cgmres.len_x*(cnt-1)), ...
                                                    u((1:cgmres.len_u)+cgmres.len_u*(cnt-1)), ...
                                                    lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cnt-1)), sys, cgmres );
    end
end