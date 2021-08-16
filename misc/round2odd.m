function [ a ] = round2odd( a )

a = round( (a - 1 ) / 2 ) * 2 + 1;

% round( (a - 1 ) / 2 ) * 2 + 1            % round to odd
% round( (a - 2 ) / 2 ) * 2 + 2            % round to even