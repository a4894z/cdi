function [ Xt ] = truncatearray( X, szt )

sz = 0.5 * size( X );
    
Xt = X(  ( sz( 1 ) + 1 - 0.5 * szt( 1 )) :  ( sz( 1 ) + 0.5 * szt( 1 )), ...
         ( sz( 2 ) + 1 - 0.5 * szt( 2 )) :  ( sz( 2 ) + 0.5 * szt( 2 )));
