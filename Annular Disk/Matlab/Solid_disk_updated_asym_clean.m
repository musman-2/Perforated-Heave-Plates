close all
clear all
clc
format long

eta = 10^(1.0);

pp = 3;

oo = 3;

n = 10;

infty = 1e7;

m = 10; 

% RTol = 1e-12;
RTol = 0;
ATol = 1e-12;

BB = zeros( m , 1 );
BB(1) = 1;

AA = zeros( m , m );
AA_alt_1 = zeros( m , m );
AA_alt_2_1 = zeros( m , m );
AA_alt_2_sin = zeros( m , m );
AA_alt_2_cos = zeros( m , m );
AA_alt = zeros( m , m );

aa_fun = @(i,j,k) k .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) .* besselj( 2 * i + 1 / 2 , k ) .* besselj( 2 * j + 1 / 2 , k );
aa_fun_r = @(i,j,k) k .* ( 1 - sqrt( ( 1 + sqrt( 1 + ( eta ./ k ).^4 ) ) ./ ( 2 * ( 1 + ( eta ./ k ) .^4 ) ) ) ) .* besselj( 2 * i + 1 / 2 , k ) .* besselj( 2 * j + 1 / 2 , k );
aa_fun_i = @(i,j,k) k .* sqrt( ( sqrt( 1 + ( eta ./ k ).^4 ) - 1 ) ./ ( 2 * ( 1 + ( eta ./ k ) .^4 ) ) ) .* besselj( 2 * i + 1 / 2 , k ) .* besselj( 2 * j + 1 / 2 , k );


aa_fun_alt_1 = @(i,j,k) ( 1 / pi ) .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) .* ( R1( 2 * i , k ) .* R1( 2 * j , k ) + R2( 2 * i - 1 , k ) .* R2( 2 * j - 1 , k ) );
aa_fun_alt_1_r = @(i,j,k) ( 1 / pi ) .* ( 1 - sqrt( ( 1 + sqrt( 1 + ( eta ./ k ).^4 ) ) ./ ( 2 * ( 1 + ( eta ./ k ) .^4 ) ) ) ) .* ( R1( 2 * i , k ) .* R1( 2 * j , k ) + R2( 2 * i - 1 , k ) .* R2( 2 * j - 1 , k ) );
aa_fun_alt_1_i = @(i,j,k) ( 1 / pi ) .* sqrt( ( sqrt( 1 + ( eta ./ k ).^4 ) - 1 ) ./ ( 2 * ( 1 + ( eta ./ k ) .^4 ) ) ) .* ( R1( 2 * i , k ) .* R1( 2 * j , k ) + R2( 2 * i - 1 , k ) .* R2( 2 * j - 1 , k ) );


aa_fun_alt_sin = @(i,j,k) - ( 1 / pi ) .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) .* ( R1( 2 * i , k ) .* R2( 2 * j - 1 , k ) + R1( 2 * j , k ) .* R2( 2 * i - 1 , k ) ) .* sin( 2 * k );
aa_fun_alt_cos = @(i,j,k) - ( 1 / pi ) .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) .* ( R1( 2 * i , k ) .* R1( 2 * j , k ) - R2( 2 * i - 1 , k ) .* R2( 2 * j - 1 , k ) ) .* cos( 2 * k );


for ii = 0 : m - 1    

    AA( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun(i,j,k) , 0 , inf , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );
    % AA( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun(i,j,k) , 0 , pi * max( oo * 10^pp * max( i + 1 , j + 1 ) , 10^pp * floor(eta) ) , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );
    % AA( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun(i,j,k) , 0 , infty , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );
    
    AA_alt_1( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun(i,j,k) , 0 , oo * pi * max( i + 1 , j + 1 ) , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );
   
    AA_alt_2_1( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun_alt_1(i,j,k) , oo * pi * max( i + 1 , j + 1 ) , inf , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );
    % AA_alt_2_1( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun_alt_1(i,j,k) , oo * pi * max( i + 1 , j + 1 ) , pi * max( oo * 10^pp * max( i + 1 , j + 1 ) , 10^pp * floor(eta) ) , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );    
    % AA_alt_2_1( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) integral( @(k) aa_fun_alt_1(i,j,k) , oo * pi * max( i + 1 , j + 1 ) , infty , 'RelTol' , RTol , 'AbsTol' , ATol ) , ii * ones( 1 , m - ii ) , ii : m - 1 );


    AA_alt_2_sin( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) int_sin( i , j , eta , pp , RTol , ATol , n , oo ) , ii * ones( 1 , m - ii ) , ii : m - 1 );

    AA_alt_2_cos( ii + 1 , ii + 1 : m ) = arrayfun( @(i,j) int_cos( i , j , eta , pp , RTol , ATol , n , oo ) , ii * ones( 1 , m - ii ) , ii : m - 1 );

end

AA = AA + transpose(AA) - diag( diag(AA) );
AA_alt_2_sin( 1 , 1 ) = 0;
AA_alt = AA_alt_1 + AA_alt_2_1 + AA_alt_2_sin + AA_alt_2_cos;
AA_alt = AA_alt + transpose(AA_alt) - diag( diag(AA_alt) );

aa = linsolve( AA , BB );
aa_alt = linsolve( AA_alt , BB );

% F = - 1i * 8 * aa(1);
% F_alt = - 1i* 8 * aa_alt(1);

aa1 = aa(1)

aa1_alt = aa_alt(1)

% F_inf_theory = 8 / 3;
% 
% F_corr_high = [ - imag(aa(1)) * eta ( real( aa(1) ) - F_inf_theory ) * eta ]
% 
% F_alt_corr_high = [ -imag( aa_alt(1) ) * eta ( real( aa_alt(1) ) - F_inf_theory ) * eta ]
% 
% F_corr_high_fit_log = [ 0.8179065673764463 + log(eta) * 0.1981125557583434 1.0154472141585995 + log(eta) * 0.2186193034391346 ]
% 
% F_corr_high_fit_Howard = 2 * 1.32 * exp( 1i * 0.27 *pi )

% F_0_theory = - 16 / eta^2;
% 
% F_corr_low = [ ( real(F) - F_0_theory) * eta imag( F ) * eta  ]
% 
% F_alt_corr_low = [ ( real( F_alt ) - F_0_theory ) * eta imag( F_alt ) * eta ]
% 
% F_corr_low_theory =  - 16 * 8 * exp( 1i * 0.25 *pi ) / ( 3 * pi )


% aa1_inf_theory = 1 / 3;
% 
% aa1_corr_high = [ ( real( aa(1) ) - aa1_inf_theory ) * eta imag(aa(1)) * eta  ]
% 
% aa1_alt_corr_high = [ ( real( aa_alt(1) ) - aa1_inf_theory ) * eta imag( aa_alt(1) ) * eta  ]
% 
% a_corr_high_fit_log = [ 1.048703 + log(eta) * 0.2109335 -1.0186035 - log(eta) * 0.1545515 ]
% 
% 
% aa1_0_theory = - 2 / eta^2;
% 
% aa1_corr_low = [ real( aa(1) )  * eta ( imag( aa(1) ) - aa1_0_theory ) * eta  ]
% 
% aa1_alt_corr_low = [ real( aa_alt(1) ) * eta ( imag( aa_alt(1) ) - aa1_0_theory ) * eta ]
% 
% aa1_corr_low_theory =  - 1i * 16 * exp( 1i * 0.25 *pi ) / ( 3 * pi )


p_o_jj = arrayfun( @(jj) sum( gamma( 1 : jj ) .* transpose( aa( 1 : jj ) ) ./ gamma( ( 1 : jj ) - 1 / 2  ) ) , 1 : m );
p_o_sh_kk_jj = zeros( floor( m / 2 - 1 / 2 ) , m );
p_o_sh_kk_jj( 1 , 1 : m - 2 ) = arrayfun( @(jj) p_o_jj( jj + 2 ) - ( p_o_jj( jj + 2 ) - p_o_jj( jj + 1 ) )^2 / ( p_o_jj( jj + 2 ) - 2 * p_o_jj( jj + 1 ) + p_o_jj( jj + 2 ) ) , 1 : m - 2 );

for kk = 2 : floor( m / 2 - 1 /2 )
    for jj = 1 : m - 2 * kk
        p_o_sh_kk_jj( kk , jj ) = p_o_sh_kk_jj( kk - 1, jj + 2 ) - ( p_o_sh_kk_jj( kk - 1, jj + 2 ) - p_o_sh_kk_jj( kk - 1, jj + 1 ) )^2 / ( p_o_sh_kk_jj( kk - 1, jj + 2 ) - 2 * p_o_sh_kk_jj( kk - 1, jj + 1 ) + p_o_sh_kk_jj( kk - 1, jj + 2 ) );
    end
end

% p_o = 1i * ( 2 / sqrt(pi) ) * sum( gamma( 1 : m ) .* transpose( aa ) ./ gamma( ( 1 : m ) - 1 / 2 ) )
p_o_sh = 1i * ( 2 / sqrt(pi) ) * p_o_sh_kk_jj( end , rem( m - 1 , 2 ) + 1 )


p_o_alt_jj = arrayfun( @(jj) sum( gamma( 1 : jj ) .* transpose( aa_alt( 1 : jj ) ) ./ gamma( ( 1 : jj ) - 1 / 2  ) ) , 1 : m );
p_o_alt_sh_kk_jj = zeros( floor( m / 2 - 1 / 2 ) , m );
p_o_alt_sh_kk_jj( 1 , 1 : m - 2 ) = arrayfun( @(jj) p_o_alt_jj( jj + 2 ) - ( p_o_alt_jj( jj + 2 ) - p_o_alt_jj( jj + 1 ) )^2 / ( p_o_alt_jj( jj + 2 ) - 2 * p_o_alt_jj( jj + 1 ) + p_o_alt_jj( jj + 2 ) ) , 1 : m - 2 );

for kk = 2 : floor( m / 2 - 1 /2 )
    for jj = 1 : m - 2 * kk
        p_o_alt_sh_kk_jj( kk , jj ) = p_o_alt_sh_kk_jj( kk - 1, jj + 2 ) - ( p_o_alt_sh_kk_jj( kk - 1, jj + 2 ) - p_o_alt_sh_kk_jj( kk - 1, jj + 1 ) )^2 / ( p_o_alt_sh_kk_jj( kk - 1, jj + 2 ) - 2 * p_o_alt_sh_kk_jj( kk - 1, jj + 1 ) + p_o_alt_sh_kk_jj( kk - 1, jj + 2 ) );
    end
end

% p_o_alt = 1i * ( 2 / sqrt(pi) ) * sum( gamma( 1 : m ) .* transpose( aa_alt ) ./ gamma( ( 1 : m ) - 1 / 2 ) )
p_o_alt_sh = 1i * ( 2 / sqrt(pi) ) * p_o_alt_sh_kk_jj( end , rem( m - 1 , 2 ) + 1 )



% p_o_inf_theory = 1i * 2 / pi
% 
% p_o_corr_high_fit_log = [ 0.9832445 + log(eta) * 0.1260717 1.1146177 + log(eta) * 0.1434451 ]
% 
% p_o_corr_high = [ real( p_o ) * eta ( imag( p_o ) - imag( p_o_inf_theory ) ) * eta ]
% p_o_sh_corr_high = [ real( p_o_sh ) * eta ( imag( p_o_sh ) - imag( p_o_inf_theory ) ) * eta ]
% 
% p_o_alt_corr_high = [ real( p_o_alt ) * eta ( imag( p_o_alt ) - imag( p_o_inf_theory ) ) * eta ]
% p_o_alt_sh_corr_high = [ real( p_o_alt_sh ) * eta ( imag( p_o_alt_sh ) - imag( p_o_inf_theory ) ) * eta ]
% 
% 
% p_o_0_theory =  ( 4 / pi ) / eta^2;
% 
% p_o_corr_low = [ ( real( p_o ) - p_o_0_theory ) * eta imag( p_o ) * eta ]
% p_o_sh_corr_low = [ ( real( p_o_sh ) - p_o_0_theory ) * eta imag( p_o_sh ) * eta ]
% 
% p_o_alt_corr_low = [ ( real( p_o_alt ) - p_o_0_theory ) * eta imag( p_o_alt ) * eta ]
% p_o_alt_sh_corr_low = [ ( real( p_o_alt_sh ) - p_o_0_theory ) * eta imag( p_o_alt_sh ) * eta ]
% 
% p_o_corr_low_theory = 32 * exp( 1i * 0.25 *pi ) / ( 3 * pi^2 )




function R_1 = R1( m , x )
    R_1 = 0;
    for i = 0 : floor( m / 2 )
        R_1 = R_1 + (-1)^i * factorial( m - i ) * gamma( 1 / 2 + m - i ) * ( x / 2 ).^( 2 * i - m) / ( factorial(i) * factorial( m - 2 * i ) * gamma( 1 / 2 + i ) );
    end
end

function R_2 = R2( m , x )
    if m == -1
        R_2 = 0;
    else
        R_2 = 0;
        for i = 0 : floor( m / 2 )
            R_2 = R_2 + (-1)^i * factorial( m - i ) * gamma( 3 / 2 + m - i ) * ( x / 2 ).^( 2 * i - m ) / ( factorial(i) * factorial( m - 2 * i ) * gamma( 3 / 2 + i ) ) ;
        end                
    end
end

function int = int_sin( i , j , eta , pp , RTol , ATol , n , oo )
        
    aa_fun_alt_sin = @(i,j,k) - ( 1 / pi ) .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) .* ( R1( 2 * i , k ) .* R2( 2 * j - 1 , k ) + R1( 2 * j , k ) .* R2( 2 * i - 1 , k ) ) .* sin( 2 * k );
    aa_fun_alt_sin_r = @(i,j,k) - ( 1 / pi ) .* ( 1 - sqrt( ( 1 + sqrt( 1 + ( eta ./ k ).^4 ) ) ./ ( 2 * ( 1 + ( eta ./ k ) .^4 ) ) ) ) .* ( R1( 2 * i , k ) .* R2( 2 * j - 1 , k ) + R1( 2 * j , k ) .* R2( 2 * i - 1 , k ) ) .* sin( 2 * k );
    aa_fun_alt_sin_i = @(i,j,k) - ( 1 / pi ) .* sqrt( ( sqrt( 1 + ( eta ./ k ).^4 ) - 1 ) ./ ( 2 * ( 1 + ( eta ./ k ) .^4 ) ) ) .* ( R1( 2 * i , k ) .* R2( 2 * j - 1 , k ) + R1( 2 * j , k ) .* R2( 2 * i - 1 , k ) ) .* sin( 2 * k );

    n1 = oo * max( i + 1 , j + 1 );
    interval = n1 * pi + ( 0 : n + 3 ) * pi / 2;
       
    int_ii = zeros( 1 , n + 3 );
    int_ii_r = zeros( 1 , n + 3 );
    int_ii_i = zeros( 1 , n + 3 );

    psi = zeros( 1 , n + 2 );
    psi_r = zeros( 1 , n + 2 );
    psi_i = zeros( 1 , n + 2 );

    F = zeros( 1 , n + 2 );
    F_r = zeros( 1 , n + 2 );
    F_i = zeros( 1 , n + 2 );

    M = zeros( n + 2 , n + 2 );
    M_r = zeros( n + 2 , n + 2 );
    M_i = zeros( n + 2 , n + 2 );

    N = zeros( n + 2 , n + 2 );
    N_r = zeros( n + 2 , n + 2 );
    N_i = zeros( n + 2 , n + 2 );

    int_ii = arrayfun( @(ii) integral( @(k) aa_fun_alt_sin(i,j,k) , interval(ii) , interval(ii + 1) , 'RelTol' , RTol , 'AbsTol' , ATol ) , 1 : n + 3 );
    int_ii_r = arrayfun( @(ii) integral( @(k) aa_fun_alt_sin_r(i,j,k) , interval(ii) , interval(ii + 1) , 'RelTol' , RTol , 'AbsTol' , ATol ) , 1 : n + 3 );
    int_ii_i = arrayfun( @(ii) integral( @(k) aa_fun_alt_sin_i(i,j,k) , interval(ii) , interval(ii + 1) , 'RelTol' , RTol , 'AbsTol' , ATol ) , 1 : n + 3 );

    F = arrayfun( @(jj) sum( int_ii( 1 : jj ) ) , 1 : n + 2 );
    F_r = arrayfun( @(jj) sum( int_ii_r( 1 : jj ) ) , 1 : n + 2 );
    F_i = arrayfun( @(jj) sum( int_ii_i( 1 : jj ) ) , 1 : n + 2 );

    M( 1 , 1 : n + 2 ) = F( 1 : n + 2 ) ./ int_ii( 2 : n + 3 ); 
    M_r( 1 , 1 : n + 2 ) = F_r( 1 : n + 2 ) ./ int_ii_r( 2 : n + 3 ); 
    M_i( 1 , 1 : n + 2 ) = F_i( 1 : n + 2 ) ./ int_ii_i( 2 : n + 3 ); 

    N( 1 , 1 : n + 2 ) = 1 ./ int_ii( 2 : n + 3 );
    N_r( 1 , 1 : n + 2 ) = 1 ./ int_ii_r( 2 : n + 3 );
    N_i( 1 , 1 : n + 2 ) = 1 ./ int_ii_i( 2 : n + 3 );
    
    for kk = 2 : n + 2
        for jj = 1 : n + 2 - kk + 1
            M( kk , jj ) = ( M( kk - 1 , jj ) - M( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            M_r( kk , jj ) = ( M_r( kk - 1 , jj ) - M_r( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            M_i( kk , jj ) = ( M_i( kk - 1 , jj ) - M_i( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            N( kk , jj ) = ( N( kk - 1 , jj ) - N( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            N_r( kk , jj ) = ( N_r( kk - 1 , jj ) - N_r( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            N_i( kk , jj ) = ( N_i( kk - 1 , jj ) - N_i( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
        end
    end
    
    int = M( n + 2 , 1 ) / N( n + 2 , 1 );
    int_r = M_r( n + 2 , 1 ) / N_r( n + 2 , 1 );
    int_i = M_i( n + 2 , 1 ) / N_i( n + 2 , 1 );
    % int = int_r + 1i * int_i;

end

function int = int_cos( i , j , eta , pp , RTol , ATol , n , oo )

    aa_fun_alt_cos = @(i,j,k) - ( 1 / pi ) .* ( 1 - k ./ sqrt( k.^2 + 1i * eta^2 ) ) .* ( R1( 2 * i , k ) .* R1( 2 * j , k ) - R2( 2 * i - 1 , k ) .* R2( 2 * j - 1 , k ) ) .* cos( 2 * k );
    aa_fun_alt_cos_r = @(i,j,k) - ( 1 / pi ) .* ( 1 - sqrt( ( 1 + sqrt( 1 + ( eta ./ k ).^4 ) ) ./ ( 2 * ( 1 + ( eta ./ k ) .^4 ) ) ) ) .* ( R1( 2 * i , k ) .* R1( 2 * j , k ) - R2( 2 * i - 1 , k ) .* R2( 2 * j - 1 , k ) ) .* cos( 2 * k );
    aa_fun_alt_cos_i = @(i,j,k) - ( 1 / pi ) .* sqrt( ( sqrt( 1 + ( eta ./ k ).^4 ) - 1 ) ./ ( 2 * ( 1 + ( eta ./ k ) .^4 ) ) ) .* ( R1( 2 * i , k ) .* R1( 2 * j , k ) - R2( 2 * i - 1 , k ) .* R2( 2 * j - 1 , k ) ) .* cos( 2 * k );

    n1 = oo * max( i + 1 , j + 1 );
    interval = n1 * pi + ( 0 : n + 3 ) * pi / 2;

    int_ii = zeros( 1 , n + 3 );
    int_ii_r = zeros( 1 , n + 3 );
    int_ii_i = zeros( 1 , n + 3 );

    psi = zeros( 1 , n + 2 );
    psi_r = zeros( 1 , n + 2 );
    psi_i = zeros( 1 , n + 2 );

    F = zeros( 1 , n + 2 );
    F_r = zeros( 1 , n + 2 );
    F_i = zeros( 1 , n + 2 );

    M = zeros( n + 2 , n + 2 );
    M_r = zeros( n + 2 , n + 2 );
    M_i = zeros( n + 2 , n + 2 );

    N = zeros( n + 2 , n + 2 );
    N_r = zeros( n + 2 , n + 2 );
    N_i = zeros( n + 2 , n + 2 );

    int_ii = arrayfun( @(ii) integral( @(k) aa_fun_alt_cos(i,j,k) , interval(ii) , interval(ii + 1) , 'RelTol' , RTol , 'AbsTol' , ATol ) , 1 : n + 3 );
    int_ii_r = arrayfun( @(ii) integral( @(k)  aa_fun_alt_cos_r(i,j,k) , interval(ii) , interval(ii + 1) , 'RelTol' , RTol , 'AbsTol' , ATol ) , 1 : n + 3 );
    int_ii_i = arrayfun( @(ii) integral( @(k) aa_fun_alt_cos_i(i,j,k) , interval(ii) , interval(ii + 1) , 'RelTol' , RTol , 'AbsTol' , ATol ) , 1 : n + 3 );

    F = arrayfun( @(jj) sum( int_ii( 1 : jj ) ) , 1 : n + 2 );
    F_r = arrayfun( @(jj) sum( int_ii_r( 1 : jj ) ) , 1 : n + 2 );
    F_i = arrayfun( @(jj) sum( int_ii_i( 1 : jj ) ) , 1 : n + 2 );

    M( 1 , 1 : n + 2 ) = F( 1 : n + 2 ) ./ int_ii( 2 : n + 3 ); 
    M_r( 1 , 1 : n + 2 ) = F_r( 1 : n + 2 ) ./ int_ii_r( 2 : n + 3 ); 
    M_i( 1 , 1 : n + 2 ) = F_i( 1 : n + 2 ) ./ int_ii_i( 2 : n + 3 ); 

    N( 1 , 1 : n + 2 ) = 1 ./ int_ii( 2 : n + 3 );
    N_r( 1 , 1 : n + 2 ) = 1 ./ int_ii_r( 2 : n + 3 );
    N_i( 1 , 1 : n + 2 ) = 1 ./ int_ii_i( 2 : n + 3 );
    
    for kk = 2 : n + 2
        for jj = 1 : n + 2 - kk + 1
            M( kk , jj ) = ( M( kk - 1 , jj ) - M( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            M_r( kk , jj ) = ( M_r( kk - 1 , jj ) - M_r( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            M_i( kk , jj ) = ( M_i( kk - 1 , jj ) - M_i( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            N( kk , jj ) = ( N( kk - 1 , jj ) - N( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            N_r( kk , jj ) = ( N_r( kk - 1 , jj ) - N_r( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
            N_i( kk , jj ) = ( N_i( kk - 1 , jj ) - N_i( kk - 1 , jj + 1 ) ) / ( 1 / interval( jj + 1 ) - 1 / interval( jj + kk + 1 ) );
        end
    end
    
    int = M( n + 2 , 1 ) / N( n + 2 , 1 );
    int_r = M_r( n + 2 , 1 ) / N_r( n + 2 , 1 );
    int_i = M_i( n + 2 , 1 ) / N_i( n + 2 , 1 );
    % int = int_r + 1i * int_i;

end