close all;
clear;

dSamplingPeriod = 0.001;

a = -1;
b = -1;

umax = 1;

tsim0 = 0;
tsim = 20;
ht = 0.001;
tf = 1.0;
alpha = 0.5;
zeta = 1000.0;
dv = 5;

x0 = [2;0];
u0 = [0.01;0.9;0.03];

len.x = length( x0 );
len.u = length( u0 );
len.lmd = len.x;

q = [ 1;10 ];
r = [ 1;0.01 ];
sf = [ 1;10 ];

tic;
sim( 'CGMRES' );
toc