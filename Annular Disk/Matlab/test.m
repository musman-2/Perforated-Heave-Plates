% close all
figHandles = findall(groot, 'Type', 'figure'); 
arrayfun(@(fh) clf(fh), figHandles);

% figure (1)
% loglog( eta , abs( abs( real(sin_integral) ) - pi / 2 ) .* eta , 'k' , 'LineWidth' , 2 )
% hold on
% loglog( eta , ( sqrt(2) / 4 ) , 'or' , 'LineWidth' , 2 )
% hold off
% 
% figure (2)
% loglog( eta , abs( abs( imag(sin_integral) ) .* eta - sqrt(2) / 4 ) , 'k' , 'LineWidth' , 2 )
% hold on
% loglog( eta , ( sqrt(3) / 20 ) ./ eta.^2 , 'or' , 'LineWidth' , 2 )
% hold off

% figure (3)
% % loglog( eta , abs( abs( real(cos_integral) ) .* eta - sqrt(2) / 8 ) , 'k' , 'LineWidth' , 2 )
% loglog( eta , 0.1302- abs( real(cos_integral) ) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( sqrt(2) / 8 ) ./ eta.^2 , 'or' , 'LineWidth' , 2 )
% hold off
% 
% figure (4)
% % loglog( eta , abs( abs( imag(cos_integral) ) .* eta - sqrt(2) / 8 ) , 'k' , 'LineWidth' , 2 )
% loglog( eta , abs( imag(cos_integral) ) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( sqrt(2) / 8 ) ./ eta.^2 , 'or' , 'LineWidth' , 2 )
% loglog( eta , 0.02 ./ eta , 'or' , 'LineWidth' , 2 )
% hold off


AA_11_r = real( AA_all( : , 1 , 1 ) );
AA_11_i = imag( AA_all( : , 1 , 1 ) );

AA_12_r = real( AA_all( : , 1 , 2 ) );
AA_12_i = imag( AA_all( : , 1 , 2 ) );

AA_22_r = real( AA_all( : , 2 , 2 ) );
AA_22_i = imag( AA_all( : , 2 , 2 ) );

AA_23_r = real( AA_all( : , 2 , 3 ) );
AA_23_i = imag( AA_all( : , 2 , 3 ) );

AA_33_r = real( AA_all( : , 3 , 3 ) );
AA_33_i = imag( AA_all( : , 3 , 3 ) );

AA_alt_11_r = real( AA_alt_all( : , 1 , 1 ) );
AA_alt_11_i = imag( AA_alt_all( : , 1 , 1 ) );

AA_alt_12_r = real( AA_alt_all( : , 1 , 2 ) );
AA_alt_12_i = imag( AA_alt_all( : , 1 , 2 ) );

AA_alt_22_r = real( AA_alt_all( : , 2 , 2 ) );
AA_alt_22_i = imag( AA_alt_all( : , 2 , 2 ) );

AA_alt_23_r = real( AA_alt_all( : , 2 , 3 ) );
AA_alt_23_i = imag( AA_alt_all( : , 2 , 3 ) );

AA_alt_33_r = real( AA_alt_all( : , 3 , 3 ) );
AA_alt_33_i = imag( AA_alt_all( : , 3 , 3 ) );

% figure (111)
% loglog( eta , AA_11_r , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.3 ) * eta.^3 , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( 1 / ( pi * sqrt(2) ) ) * eta , 'or' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-10 1e2])
% % axis([1e1 1e4])
% 
% figure (112)
% loglog( eta , AA_11_i , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.5 ) * eta.^2 , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( 1 / ( pi * sqrt(2) ) ) * eta , 'or' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-7 1e2])
% % axis([1e1 1e4])
% 
% figure (113)
% % loglog( eta , 0.3 - AA_11_r' ./ eta.^3, 'k' , 'LineWidth' , 2 )
% loglog( eta , abs( AA_alt_11_r' - eta / ( pi * sqrt(2) ) ) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.2 ) * eta , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( 0.056 ) ./ eta , 'ob' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-10 1e2])
% 
% figure (114)
% % loglog( eta , 0.5 - AA_11_i' ./ eta.^2 , 'k' , 'LineWidth' , 2 )
% loglog( eta , abs( AA_alt_11_i' - eta / ( pi * sqrt(2) ) ) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.3 ) * eta , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( 0.056 ) ./ eta , 'ob' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-7 1e2])
% 
% 
% figure (121)
% loglog( eta , abs(AA_12_r) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.3 ) * eta.^4 , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( 1 / ( pi * sqrt(2) ) ) * eta , 'or' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-10 1e2])
% 
% figure (122)
% loglog( eta , abs(AA_12_i) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.25 ) * eta.^5 , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( 1 / ( pi * sqrt(2) ) ) * eta , 'or' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-7 1e2])

a_AA_alt_12_r = ( ( 3 / 2 - abs( abs(AA_alt_12_r(4)) - eta(4) / ( pi * sqrt(2) ) ) ) * eta(4) - ( 3 / 2 - abs( abs(AA_alt_12_r(3)) - eta(3) / ( pi * sqrt(2) ) ) ) * eta(3) ) / ( log( eta(4) ) - log( eta(3) ) );
b_AA_alt_12_r = ( 3 / 2 - abs( abs(AA_alt_12_r(4) ) - eta(4) / ( pi * sqrt(2) ) ) ) * eta(4) - a_AA_alt_12_r * log( eta(4) );

figure (123)
% loglog( eta , 0.3 - AA_12_r' ./ eta.^3, 'k' , 'LineWidth' , 2 )
loglog( eta , 3 / 2 - abs( abs( AA_alt_12_r' ) - eta / ( pi * sqrt(2) ) ) , 'k' , 'LineWidth' , 2 )
hold on
% loglog( eta , ( 0.225 ) .* eta , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( 0.056 ) .* eta , 'ob' , 'LineWidth' , 2 )
loglog( eta , ( b_AA_alt_12_r + a_AA_alt_12_r * log(eta) ) ./ eta , 'ob' , 'LineWidth' , 2 )
hold off
% axis([1e-3 1e2 1e-10 1e2])

a_AA_alt_12_i = ( abs( abs(AA_alt_12_i(4)) - eta(4) / ( pi * sqrt(2) ) ) * eta(4) - abs( abs(AA_alt_12_i(3)) - eta(3) / ( pi * sqrt(2) ) ) * eta(3) ) / ( log( eta(4) ) - log( eta(3) ) );
b_AA_alt_12_i = abs( abs(AA_alt_12_i(4)) - eta(4) / ( pi * sqrt(2) ) ) * eta(4) - a_AA_alt_12_i * log( eta(4) );


figure (124)
% loglog( eta , 0.5 - AA_12_i' ./ eta.^2 , 'k' , 'LineWidth' , 2 )
loglog( eta , abs( abs( AA_alt_12_i' ) - eta / ( pi * sqrt(2) )  ) , 'k' , 'LineWidth' , 2 )
hold on
% loglog( eta , ( 0.3 ) * eta , 'ob' , 'LineWidth' , 2 )
loglog( eta , ( b_AA_alt_12_i + a_AA_alt_12_i * log(eta) ) ./ eta , 'ob' , 'LineWidth' , 2 )
hold off
% axis([1e-3 1e2 1e-7 1e2])


% figure (221)
% loglog( eta , abs(AA_22_r) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.3 ) * eta.^4 , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( 1 / ( pi * sqrt(2) ) ) * eta , 'or' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-10 1e2])
% 
% figure (222)
% loglog( eta , abs(AA_22_i) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.25 ) * eta.^5 , 'ob' , 'LineWidth' , 2 )
% loglog( eta , (  1 / ( pi * sqrt(2) ) ) * eta , 'or' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-7 1e2])
% 
% 
% a_AA_alt_22_r = ( abs( abs(AA_alt_22_r(4)) - eta(4) / ( pi * sqrt(2) ) ) * eta(4) - abs( abs(AA_alt_22_r(3)) - eta(3) / ( pi * sqrt(2) ) ) * eta(3) ) / ( log( eta(4) ) - log( eta(3) ) );
% b_AA_alt_22_r = abs( abs(AA_alt_22_r(4)) - eta(4) / ( pi * sqrt(2) ) ) * eta(4) - a_AA_alt_22_r * log( eta(4) );
% 
% a_AA_alt_22_i = ( abs( abs(AA_alt_22_i(4)) - eta(4) / ( pi * sqrt(2) ) ) * eta(4) - abs( abs(AA_alt_22_i(3)) - eta(3) / ( pi * sqrt(2) ) ) * eta(3) ) / ( log( eta(4) ) - log( eta(3) ) );
% b_AA_alt_22_i = abs( abs(AA_alt_22_i(4)) - eta(4) / ( pi * sqrt(2) ) ) * eta(4) - a_AA_alt_22_i * log( eta(4) );
% 
% figure (223)
% % loglog( eta , 0.3 - AA_12_r' ./ eta.^3, 'k' , 'LineWidth' , 2 )
% loglog( eta , abs( abs(AA_alt_22_r') - eta / ( pi * sqrt(2) ) ) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 1 / ( pi * sqrt(2) ) ) * eta , 'ob' , 'LineWidth' , 2 )
% loglog( eta ,  1 , 'ob' , 'LineWidth' , 2 )
% loglog( eta ,  ( b_AA_alt_22_r + a_AA_alt_22_r * log(eta) ) ./ eta , 'or' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-10 1e2])
% 
% figure (224)
% % loglog( eta , 0.5 - AA_12_i' ./ eta.^2 , 'k' , 'LineWidth' , 2 )
% loglog( eta , abs( abs(AA_alt_22_i') - eta / ( pi * sqrt(2) ) ) , 'k' , 'LineWidth' , 2 )
% hold on
% % loglog( eta , ( 0.3 ) * eta , 'ob' , 'LineWidth' , 2 )
% loglog( eta , ( b_AA_alt_22_i + a_AA_alt_22_i * log(eta) ) ./ eta , 'ob' , 'LineWidth' , 2 )
% hold off
% % axis([1e-3 1e2 1e-7 1e2])