function Phix = dPhidx( x, sf )
    Phix = [ x(1) * sf(1);x(2) * sf(2) ];
end