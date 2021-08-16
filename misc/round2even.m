function [ a ] = round2even( a )

a = round( (a - 2 ) / 2 ) * 2 + 2;

% round( (a - 1 ) / 2 ) * 2 + 1            % round to odd
% round( (a - 2 ) / 2 ) * 2 + 2            % round to even