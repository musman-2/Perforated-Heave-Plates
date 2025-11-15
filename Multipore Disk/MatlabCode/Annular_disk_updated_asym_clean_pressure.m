%%

clear all
clc
format long

eta = 10^(-2.0);

pp = 3;

oo = 3;

infty_1 = 1e7;
% infty_2 = 1e5;

m = 10; 

RTol = 1e-12;
% ATol = 1e-15;

BB = zeros( m , 1 );
BB(1) = 1;

AA = zeros( m , m );

aa_fun = @( i , j , k ) besselj( 2 * i + 1 / 2 , k ) .* besselj( 2 * j + 1 / 2 , k ) ./ ( k .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) );
aa_fun_p = @( i , rho , k ) besselj( 2 * i + 1 / 2 , k ) .* besselj( 0 , k * rho ) ./ ( sqrt( k ) .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) );

for ii = 0 : m - 1    

    AA( ii + 1 , ii + 1 : m ) = arrayfun( @( i , j ) integral( @(k) aa_fun(i,j,k) , 0 , infty_1 , 'RelTol' , RTol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );

end

AA = AA + transpose(AA) - diag( diag(AA) );
aa = linsolve( AA , BB );

rho = 10.^( 0.125 : 0.125 : 2 );

infty_2 = 3e4 ./ rho.^1.0;
infty_2( 1 : 10 ) = 5e4 ./ rho( 1 : 10 ).^1.05;

for ii = 1 : length( rho )

    int_p = arrayfun( @( i ) integral( @(k) aa_fun_p( i , rho( ii ) , k ) , 0 , infty_2( ii ) , 'RelTol' , RTol ) , 0 : m - 1 );

    p_o_jj = arrayfun( @( kk ) sum( transpose( aa( 1 : kk ) ) .* int_p( 1 : kk ) ) , 1 : m );

    p_o_sh_kk_jj = zeros( floor( m / 2 - 1 / 2 ) , m );
    p_o_sh_kk_jj( 1 , 1 : m - 2 ) = arrayfun( @(jj) p_o_jj( jj + 2 ) - ( p_o_jj( jj + 2 ) - p_o_jj( jj + 1 ) )^2 / ( p_o_jj( jj + 2 ) - 2 * p_o_jj( jj + 1 ) + p_o_jj( jj + 2 ) ) , 1 : m - 2 );
    
    for kk = 2 : floor( m / 2 - 1 /2 )

        for jj = 1 : m - 2 * kk
            p_o_sh_kk_jj( kk , jj ) = p_o_sh_kk_jj( kk - 1, jj + 2 ) - ( p_o_sh_kk_jj( kk - 1, jj + 2 ) - p_o_sh_kk_jj( kk - 1, jj + 1 ) )^2 / ( p_o_sh_kk_jj( kk - 1, jj + 2 ) - 2 * p_o_sh_kk_jj( kk - 1, jj + 1 ) + p_o_sh_kk_jj( kk - 1, jj + 2 ) );
        end

    end
    
    p_o_sh(ii) = sqrt( 2 / pi ) * p_o_sh_kk_jj( end , rem( m - 1 , 2 ) + 1 );

    % p_o_sh(ii) = sqrt( 2 / pi ) * p_o_jj( m );

end


clipboard( 'copy' , sprintf( '%.8f\t%.8f\n' , [real(p_o_sh(:)) imag(p_o_sh(:))]' ) );

%%

close all

figure(1)

subplot(3,2,1)
loglog( rho , abs( real( p_o_sh ) ) , 'k' , 'LineWidth', 1 )
hold on
loglog( rho , abs( ( 2 / pi ) * ( asin( 1./ rho ) - 1 ./ sqrt( rho.^2 - 1 ) ) ) , 'bo' , 'LineWidth', 1 )
loglog( rho , ( 2 / pi ) * asin( 1 ./ rho ) , 'ro' , 'LineWidth', 1 )
hold off

subplot(3,2,2)
loglog( rho , abs( imag( p_o_sh ) ) , 'k' , 'LineWidth', 1 )
hold on
loglog( rho , ( 2 / pi ) * ( asin( 1./ rho ) - 1 ./ sqrt( rho.^2 - 1 ) ) , 'bo' , 'LineWidth', 1 )
loglog( rho , ( 2 / pi ) * asin( 1 ./ rho ) , 'ro' , 'LineWidth', 1 )
hold off

subplot(3,2,3)
semilogx( rho , 1e6 * real( p_o_sh ) , 'k' , 'LineWidth', 1 )
hold on
semilogx( rho , 1e6 * ( 2 / pi ) * ( asin( 1./ rho ) - 1 ./ sqrt( rho.^2 - 1 ) ) , 'bo' , 'LineWidth', 1 )
% semilogx( rho , ( 2 / pi ) * asin( 1 ./ rho ) , 'ro' , 'LineWidth', 1 )
hold off

subplot(3,2,4)
semilogx( rho , imag( p_o_sh ) , 'k' , 'LineWidth', 1 )
hold on
semilogx( rho , ( 1 ./ sqrt( rho.^2 - 1 ) - asin( 1./ rho ) ) , 'bo' , 'LineWidth', 1 )
semilogx( rho , ( 2 / pi ) * asin( 1 ./ rho ) , 'ro' , 'LineWidth', 1 )
hold off

subplot(3,2,5)
loglog( rho , abs( real( p_o_sh ) / real( p_o_sh(1) ) ) , 'k' , 'LineWidth', 1 )
hold on
loglog( rho , ( asin( 1./ rho ) - 1 ./ sqrt( rho.^2 - 1 ) ) / ( asin( 1/ rho(1) ) - 1 ./ sqrt( rho(1)^2 - 1 ) ) , 'bo' , 'LineWidth', 1 )
loglog( rho , asin( 1 ./ rho ) / asin( 1 / rho(1) ) , 'ro' , 'LineWidth', 1 )
hold off

subplot(3,2,6)
loglog( rho , abs( imag( p_o_sh ) / imag( p_o_sh(1) ) ) , 'k' , 'LineWidth', 1 )
hold on
loglog( rho , ( asin( 1./ rho ) - 1 ./ sqrt( rho.^2 - 1 ) ) / ( asin( 1/ rho(1) ) - 1 ./ sqrt( rho(1)^2 - 1 ) ) , 'bo' , 'LineWidth', 1 )
loglog( rho , asin( 1 ./ rho ) / asin( 1 / rho(1) ) , 'ro' , 'LineWidth', 1 )
hold off

%%

% clc
% 
% r_0 = 100;
% 
% test_1 = @( r ) integral( @(k) sin( k ) .* besselj( 0 , k * r ) ./ k , 0 , inf , 'RelTol' , 0 , 'AbsTol' , 1e-12 );
% test_11 = @( r ) asin( 1 / r );
% test_2 = @( r ) integral( @(k) cos( k ) .* besselj( 0 , k * r ) , 0 , 1e2 , 'RelTol' , 1e-12 );
% test_22 = @( r ) 1 / sqrt( r^2 - 1 );
% 
% test_1( r_0 )
% test_11( r_0 )
% test_2( r_0 )
% test_22( r_0 )


