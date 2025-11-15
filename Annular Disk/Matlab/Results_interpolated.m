close all
clear all
clc
format long

load( 'alpha.mat' );
load( 'p.mat' );
load( 'beta.mat' );

N = 1000;
 
% eta = 10.^linspace( -5 , 5 , N + 1 );
% eta_eps = 10.^linspace( -5 , 5 , N + 1 );
% eps = linspace( 0 , 1 / 2 , N + 1 );

eta = ones( 36 , 1 ) * sqrt(10.^( -4 : 1 : 6 ) );
% eta = ones( length(eps) , 1 ) * 10.^( 2 : 0.5 : 6 ) ;
% eta = ones( length(eps) , 1 ) * 10.^( 0.5 : 0.5 : 1.5 ) ;
% eps = ( 0.08 : 0.001 : 0.5 )';
% eps = ( 0.02 : 0.0001 : 0.08 )';
% eps = linspace( 0 , 1 , N + 1 );
eps= [
0
10^-6
10^-5.5
10^-5
10^-4.5
10^-4
10^-3.5
10^-3
10^-2.75
10^-2.5
10^-2.25
0.01
0.02
0.03
0.04
0.06
0.08
0.1
0.15
0.2
0.25
0.3
0.35
0.4
0.45
0.5
0.55
0.6
0.65
0.7
0.75
0.8
0.85
0.9
0.95
0.99];

eps = eps * ones( 1 , size( eta , 2 ) );

eta_eps = eta .* eps;

alpha_0 = alpha0( eta , alpha_0_eta );
pp_o = po( eta , pp_o_eta );
beta_0 = beta0( eta_eps , beta_0_eta_eps );

F_0 = - 1i * 8 * eta.^2 .* alpha_0;
F_1 = - 1i * 8 * eta.^2 .* eps .* pp_o .* pp_o .* beta_0;

a_0 = 8 * real(alpha_0);
b_0 = -8 * imag(alpha_0);

a_1 = 8 * eps .* real(pp_o .* pp_o .* beta_0);
b_1 = -8 * eps .* imag(pp_o .* pp_o .* beta_0);

a = a_0 + a_1;
b = b_0 + b_1;

a_norm = 1 + eps .* real(pp_o .* pp_o .* beta_0) ./ real(alpha_0);
a_norm( 1 , : ) = 1;
b_norm = 1 + eps .* imag(pp_o .* pp_o .* beta_0) ./ imag(alpha_0);
b_norm( 1 , : ) = 1;

F_norm = sqrt( 1 + ( 2 * ( a_0 .* a_1 + b_0 .* b_1 ) + a_1.^2 + b_1.^2 ) ./ ( a_0.^2 + b_0.^2 ) );
F_norm( 1 , : ) = 1;
% F_norm_alt = sqrt( ( a.^2 + b.^2 ) ./ ( a_0.^2 + b_0.^2 ) );
% F_norm_alt( 1 , : ) = 1;

phi_norm = atan( a ./ b ) ./ atan( a_0 ./ b_0 );
phi_norm( 1 , : ) = 1;

phi = phi_norm .* atan( a_0 ./ b_0 );

% tiledlayout(1,2)
% 
% nexttile
% plot( eps , a_norm , 'k' )
% % hold on
% % plot( eps , a_norm_alt_asym , 'o' )
% % hold off
% 
% nexttile
% plot( eps , b_norm , 'k' )
% % hold on
% % plot( eps , b_norm_alt_asym , 'o' )
% % hold off


function alpha = alpha0( x , table )

alpha_real = ( x < 0.01 ) .* ( 16 ./ ( 3 * sqrt(2) * pi * x ) ) ...
      + ( x > 1000 ) .* ( 1 / 3 + ( 0.2124754 * log(x) + 1.0436181 ) ./ x ) ...
      + ( ( 0.01 <= x ) & ( x <= 1000) ) .* ( interp1( log10(table.eta) , table.alpha_0_real , log10(x) , 'spline') );

alpha_imag = ( x < 0.01 ) .*  ( -2 ./ x.^2 - 16 ./ ( 3 * sqrt(2) * pi * x ) ) ...
      + ( x > 1000 ) .* ( ( -0.1543773 * log(x) - 1.0193142 ) ./ x ) ...
      + ( ( 0.01 <= x ) & ( x <= 1000) ) .* ( interp1( log10(table.eta) , table.alpha_0_imag , log10(x) , 'spline') );

alpha = alpha_real + 1i * alpha_imag;
  
end


function p = po( x , table )

p_real = ( x < 0.01 ) .* ( 4 ./ ( pi * x.^2 ) + 32 ./ ( 3 * sqrt(2) * pi^2 * x ) ) ...
      + ( x > 1000 ) .* ( ( 0.1259719 * log(x) + 0.9836419 ) ./ x ) ...
      + ( ( 0.01 <= x ) & ( x <= 1000) ) .* ( interp1( log10(table.eta) , table.pp_o_real , log10(x) , 'spline') );

p_imag = ( x < 0.01 ) .*  ( 32 ./ ( 3 * sqrt(2) * pi^2 * x ) ) ...
      + ( x > 1000 ) .* ( 2 / pi + ( 0.1442063 * log(x) + 1.1120969 ) ./ x ) ...
      + ( ( 0.01 <= x ) & ( x <= 1000) ) .* ( interp1( log10(table.eta) , table.pp_o_imag , log10(x) , 'spline') );

p = p_real + 1i * p_imag;
  
end


function beta = beta0( x , table )

beta_real = ( x < 0.01 ) .* ( x.^4 / 20 ) ...
      + ( x > 1000 ) .* ( 1 + ( -0.2270745 * log(x) - 1.1858204 ) ./ x ) ...
      + ( ( 0.01 <= x ) & ( x <= 1000) ) .* ( interp1( log10(table.eta_eps) , table.beta_0_real , log10(x) , 'spline') );

beta_imag = ( x < 0.01 ) .*  ( x.^2 / 6 ) ...
      + ( x > 1000 ) .* ( ( 0.2369283 * log(x) + 0.7801902 ) ./ x ) ...
      + ( ( 0.01 <= x ) & ( x <= 1000) ) .* ( interp1( log10(table.eta_eps) , table.beta_0_imag , log10(x) , 'spline') );

beta = beta_real + 1i * beta_imag;
  
end