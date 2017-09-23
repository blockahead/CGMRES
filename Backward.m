%% Œ»İ‚©‚çT•b–¢—ˆ‚Ü‚Å‚Ì”º•Ï”‚Ì—\‘ªiEuler‹ß—j
function lmd = Backward( x, u, T, sys, cgmres )
    dt = T / cgmres.dv;
    
    lmd = zeros( cgmres.len_lmd * cgmres.dv, 1 );
    lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cgmres.dv-1)) = dPhidx( x((1:cgmres.len_x)+cgmres.len_x*(cgmres.dv-1)), cgmres );
    
    for cnt = cgmres.dv-1:-1:1
        lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cnt-1)) = lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cnt)) ...
                                                            + dHdx( x((1:cgmres.len_x)+cgmres.len_x*(cnt)), lmd((1:cgmres.len_lmd)+cgmres.len_lmd*(cnt)), sys, cgmres ) * dt;
    end
end