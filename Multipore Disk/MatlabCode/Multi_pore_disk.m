clc
clear
close all

load( 'alpha.mat' );
load( 'p.mat' );
p_sd = load( 'p_sd.mat' );
load( 'beta.mat' );
load( 'max.mat' );
p_hiw = load( 'p_hiw.mat' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = linspace( 0, 2 * pi , 100);

N = 6;
R = 1 / 2;
eps = 0.05;
eta = 1e2;
eta_eps = eps * eta;

% x( 1 ) = 0;     x( 2 : N ) = R * cos( ( 0 : N - 2 ) * 2 * pi / ( N - 1 ) );
% y( 1 ) = 0;     y( 2 : N ) = R * sin( ( 0 : N - 2 ) * 2 * pi / ( N - 1 ) );

x( 1 : N ) = R * cos( ( 0 : N - 1 ) * 2 * pi / N );
y( 1 : N ) = R * sin( ( 0 : N - 1 ) * 2 * pi / N );

rho( 1 : N ) = sqrt( x.^2 + y.^2 );

distance = pdist2( [ x' y' ] , [ x' y' ] ) / eps;
distance = distance - diag( diag( distance ) ) + diag( ones( N , 1 ) ) / 2;
 
F_0_star = 8 * 1i * alpha0( eta , alpha_0_eta );
a_0 = 8 * imag( F_0_star );
b_0 = 8 * real( F_0_star );


p_0_star_n = p0( eta * ones( 1 , N ) , rho , p_sd );

F_1_star_1 = 8 * 1i * eps * beta0( eta_eps , beta_0_eta_eps ) * sum( p_0_star_n.^2 );
a_1_1 = 8 * imag( F_1_star_1 );
b_1_1 = 8 * real( F_1_star_1 );

b = ones( N , 1 );

A = p1( eta_eps * ones( N , N ) , distance , p_hiw ) .* ...
repmat( p0( eta * ones( 1 , N ) , rho , p_sd ) , N , 1 ) ./ repmat( p0( eta , rho , p_sd ) , 1 , N );

A = A - diag( diag( A ) ) + diag( ones( N , 1 ) );
Q_n_norm = linsolve( A , b );
 
 
F_1_star_2 = 8 * 1i * eps * beta0( eta_eps , beta_0_eta_eps ) * sum( p_0_star_n.^2 .* Q_n_norm.' );
a_1_2 = 8 * imag( F_1_star_2 );
b_1_2 = 8 * real( F_1_star_2 );

a_norm_0 = 1;
a_norm_1 = 1 + a_1_1 / a_0;
a_norm_2 = 1 + a_1_2 / a_0;
a_norm_sh = ( a_norm_2 * a_norm_0 - a_norm_1^2 ) / ( a_norm_2 - 2 * a_norm_1 + a_norm_0 );

b_norm_0 = 1;
b_norm_1 = 1 + b_1_1 / b_0
b_norm_2 = 1 + b_1_2 / b_0
b_norm_sh = ( b_norm_2 * b_norm_0 - b_norm_1^2 ) / ( b_norm_2 - 2 * b_norm_1 + b_norm_0 )

F_0_s = - 16;

p2_0_n_s = 1 ./ ( 1 - rho.^2 );

b = ones( N , 1 );
A_s = ( 1 - ( 2 / pi ) * ( 1 ./ sqrt( distance.^2 - 1 ) + atan( sqrt( distance.^2 - 1 ) ) ) ) .* sqrt( repmat( 1 - rho'.^2 , 1 , N ) ./ repmat( 1 - rho.^2 , N , 1 ) ); 
A_s = A_s - diag( diag( A_s ) ) + diag( ones( N , 1 ) );
Q_n_norm_s = linsolve( A_s , b );

F_1_s_0 = ( 64 / ( 3 * pi^2 ) ) * eps^3 * sum( p2_0_n_s );
F_1_s_1 = ( 64 / ( 3 * pi^2 ) ) * eps^3 * sum( p2_0_n_s .* Q_n_norm_s' );

b_norm_s_0 = 1;
b_norm_s_1 = 1 + F_1_s_0 / F_0_s;
b_norm_s_2 = 1 + F_1_s_1 / F_0_s;
b_norm_s_sh = ( b_norm_s_2 * b_norm_s_0 - b_norm_s_1^2 ) / ( b_norm_s_2 - 2 * b_norm_s_1 + b_norm_s_0 );


F_0_p = - 8 / 3;

p2_0_n_p = 1 - rho.^2;

b = ones( N , 1 );
A_p = ( 2 / pi ) * asin( 1 ./ distance ) .* sqrt( repmat( 1 - rho.^2 , N , 1 ) ./ repmat( 1 - rho'.^2 , 1 , N ) ); 
A_p = A_p - diag( diag( A_p ) ) + diag( ones( N , 1 ) );
Q_n_norm_p = linsolve( A_p , b );

F_1_p_0 = ( 32 / pi^2 ) * eps * sum( p2_0_n_p );
F_1_p_1 = ( 32 / pi^2 ) * eps * sum( p2_0_n_p .* Q_n_norm_p' );

a_norm_p_0 = 1;
a_norm_p_1 = 1 + F_1_p_0 / F_0_p;
a_norm_p_2 = 1 + F_1_p_1 / F_0_p;
% mean( [ a_norm_p_1 , a_norm_p_2 ] )
% a_norm_p_sh = ( 1 + F_1_p_1 / F_0_p - ( 1 + F_1_p_0 / F_0_p )^2 ) / ( ( F_1_p_1 - 2 * F_1_p_0 ) / F_0_p )
a_norm_p_sh = ( a_norm_p_2 * a_norm_p_0 - a_norm_p_1^2 ) / ( a_norm_p_2 - 2 * a_norm_p_1 + a_norm_p_0 );
a_norm_p = [ a_norm_p_1 a_norm_p_2 a_norm_p_sh ];
% clipboard('copy', sprintf('%.4f\t%.4f\t%.4f', result));


hold on
plot( x , y , 'or' , 'LineWidth', 1 , 'MarkerSize', 800*eps )
plot( cos(theta) , sin(theta) , 'k' , 'LineWidth', 1 )
hold off
xlim([-1 1])
ylim([-1 1])
axis equal
grid on

%%

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


function psd = p0( eta , rho , table )

rho_in_interp = ( 0 : 0.05 : 0.9 );
rho_in_interp = [ rho_in_interp , 1 - 10.^(-1.25 : -0.25 : -2 ) ];
eta_interp = 10.^( -2 : 0.25 : 3 );
rho_in_interp_sd = repmat( rho_in_interp' , 1 , length( eta_interp ) );
eta_interp_sd = repmat( eta_interp , length( rho_in_interp ) , 1 );

p0_real = interp2( log10( eta_interp_sd ) , rho_in_interp_sd , table.p_sd_r , log10( eta ) , rho , 'spline' );

p0_imag = interp2( log10( eta_interp_sd ) , rho_in_interp_sd , table.p_sd_i , log10( eta ) , rho , 'spline' );

psd = p0_real + 1i * p0_imag;
  
end


function phiw = p1( eta_eps , rho , table )

rho_out_interp = 10.^( 0.125 : 0.125 : 2 );
eta_eps_interp = 10.^( -2 : 0.25 : 3 );
rho_out_interp_hiw = repmat( rho_out_interp' , 1 , length( eta_eps_interp ) );
eta_eps_interp_hiw = repmat( eta_eps_interp , length( rho_out_interp ) , 1 );

p1_real = interp2( log10( eta_eps_interp_hiw ) , log10( rho_out_interp_hiw ) , table.p_hiw_r , log10( eta_eps ) , log10( rho ) , 'spline' );

p1_imag = interp2( log10( eta_eps_interp_hiw ) , log10( rho_out_interp_hiw ) , table.p_hiw_i , log10( eta_eps ) , log10( rho ) , 'spline' );

phiw = p1_real + 1i * p1_imag;
  
end







