%% Œ»İ‚©‚çT•b–¢—ˆ‚Ü‚Å‚Ìó‘Ô‚Ì—\‘ªiEuler‹ß—j
function x = Forward( x0, u, T, dv, a, b, len )
    dt = T / dv;
    
    x = zeros( len.x * dv, 1 );
    x(1:len.x) = x0;

    for cnt = 1 : dv-1
       x((1:len.x)+len.x*(cnt)) = x((1:len.x)+len.x*(cnt-1)) + Func( x((1:len.x)+len.x*(cnt-1)), u((1:len.u)+len.u*(cnt-1)), a, b ) * dt; 
    end
end