%% Œ»İ‚©‚çT•b–¢—ˆ‚Ü‚Å‚Ìó‘Ô‚Ì—\‘ªiEuler‹ß—j
function x = Forward( x0, u, T, sys, cgmres )
    dt = T / cgmres.dv;
    
    x = zeros( cgmres.len_x * cgmres.dv, 1 );
    x(1:cgmres.len_x) = x0;

    for cnt = 1 : cgmres.dv-1
       x((1:cgmres.len_x)+cgmres.len_x*(cnt)) = x((1:cgmres.len_x)+cgmres.len_x*(cnt-1)) ...
                                                + Func( x((1:cgmres.len_x)+cgmres.len_x*(cnt-1)), u((1:cgmres.len_u)+cgmres.len_u*(cnt-1)), sys ) * dt; 
    end
end