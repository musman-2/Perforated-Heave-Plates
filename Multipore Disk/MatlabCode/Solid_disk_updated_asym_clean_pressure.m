%%

close all
clear all
clc
format long

eta = 10^(2.0);

pp = 3;

oo = 3;

infty = 1e7;

m = 10; 

% RTol = 1e-12;
RTol = 0;
ATol = 1e-12;

BB = zeros( m , 1 );
BB(1) = 1;

AA = zeros( m , m );

aa_fun = @( i , j , k ) k .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) .* besselj( 2 * i + 1 / 2 , k ) .* besselj( 2 * j + 1 / 2 , k );

for ii = 0 : m - 1    

    AA( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun(i,j,k) , 0 , inf , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );
    % AA( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun(i,j,k) , 0 , pi * max( oo * 10^pp * max( i + 1 , j + 1 ) , 10^pp * floor(eta) ) , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );
    % AA( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun(i,j,k) , 0 , infty , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );

end

AA = AA + transpose(AA) - diag( diag(AA) );
aa = linsolve( AA , BB );

rho = ( 0 : 0.05 : 0.9 );
rho = [ rho , 1 - 10.^(-1.25 : -0.25 : -2 ) ];

for ii = 1 : length( rho )

    % int_p = arrayfun( @( aa , bb ) hypergeom( [ aa , bb ] , 1 , rho(ii)^2 ) , 1 : m , 3 / 2 - ( 1 : m ) );
    int_p = arrayfun( @( aa , bb ) ( 1 - rho(ii)^2 ).^( 1 - aa - bb ) .* hypergeom( [ 1 - aa , 1 - bb ] , 1 , rho(ii)^2 ) , 1 : m , 3 / 2 - ( 1 : m ) );

    p_o_jj = arrayfun( @( kk ) sum( gamma( 1 : kk ) .* transpose( aa( 1 : kk ) ) .* int_p( 1 : kk ) ./ gamma( ( 1 : kk ) - 1 / 2  ) ) , 1 : m );
    
    p_o_sh_kk_jj = zeros( floor( m / 2 - 1 / 2 ) , m );
    p_o_sh_kk_jj( 1 , 1 : m - 2 ) = arrayfun( @(jj) p_o_jj( jj + 2 ) - ( p_o_jj( jj + 2 ) - p_o_jj( jj + 1 ) )^2 / ( p_o_jj( jj + 2 ) - 2 * p_o_jj( jj + 1 ) + p_o_jj( jj + 2 ) ) , 1 : m - 2 );
    
    for kk = 2 : floor( m / 2 - 1 /2 )

        for jj = 1 : m - 2 * kk
            p_o_sh_kk_jj( kk , jj ) = p_o_sh_kk_jj( kk - 1, jj + 2 ) - ( p_o_sh_kk_jj( kk - 1, jj + 2 ) - p_o_sh_kk_jj( kk - 1, jj + 1 ) )^2 / ( p_o_sh_kk_jj( kk - 1, jj + 2 ) - 2 * p_o_sh_kk_jj( kk - 1, jj + 1 ) + p_o_sh_kk_jj( kk - 1, jj + 2 ) );
        end

    end
    
    p_o_sh(ii) = ( 2 * 1i / sqrt( pi ) ) * p_o_sh_kk_jj( end , rem( m - 1 , 2 ) + 1 );

end


clipboard( 'copy' , sprintf( '%.8f\t%.8f\n' , [real(p_o_sh(:)) imag(p_o_sh(:))]' ) );

%%

% A_p_r = 0.126;
% B_p_r = 0.984;
% A_p_i = 0.144;
% B_p_i = 1.112;
% 
% figure(1)
% hold on
% plot( rho , eta^2 * real( p_o_sh ) , 'k' , 'LineWidth', 1 )
% % plot( rho , ( 4 / pi + eta * 16 * sqrt( 2 ) / ( 3 * pi^2 ) ) ./ sqrt( 1 - rho.^2 ) , 'ro' , 'LineWidth', 1 )
% % plot( rho , ( A_p_r * eta * log(eta) + eta * B_p_r ) ./ ( 1 - rho.^2 ).^( 5 / 16 ) , 'ro' , 'LineWidth', 1 )
% plot( rho , ( A_p_r * eta * log(eta) + eta * B_p_r ) ./ ( 1 - rho.^2 ).^( 8 / 16 ) , 'ro' , 'LineWidth', 1 )
% hold off
% 
% figure(2)
% hold on
% % plot( rho , imag( p_o_sh ) , 'k' , 'LineWidth', 1 )
% % plot( rho , ( 2 / pi) * sqrt( 1 - rho.^2 ) , 'ro' , 'LineWidth', 1 )
% plot( rho , eta^2 * imag( p_o_sh ) , 'k' , 'LineWidth', 1 )
% % plot( rho , eta * ( 16 * sqrt( 2 ) / ( 3 * pi^2 ) ) ./ sqrt( 1 - rho.^2 ) , 'ro' , 'LineWidth', 1 )
% plot( rho , ( eta^2 *2 / pi + A_p_i * eta * log(eta) + eta * B_p_i ) * sqrt( 1 - rho.^2 ) , 'ro' , 'LineWidth', 1 )
% hold off
% 
% figure(11)
% hold on
% plot( rho , real( p_o_sh ) / real( p_o_sh(1) ) , 'k' , 'LineWidth', 1 )
% % plot( rho , 1 ./ sqrt( 1 - rho.^2 ) , 'ro' , 'LineWidth', 1 )
% plot( rho , 1 ./ ( 1 - rho.^2 ).^( 5 / 16 ) , 'ro' , 'LineWidth', 1 )
% hold off
% 
% figure(22)
% hold on
% plot( rho , imag( p_o_sh ) / imag( p_o_sh(1) ) , 'k' , 'LineWidth', 1 )
% plot( rho , sqrt( 1 - rho.^2 ) , 'ro' , 'LineWidth', 1 )
% % plot( rho , 1 ./ sqrt( 1 - rho.^2 ) , 'ro' , 'LineWidth', 1 )
% hold off
