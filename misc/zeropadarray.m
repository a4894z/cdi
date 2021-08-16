function [  X ] = zeropadarray( X, szp )

% zeropad rows:
if szp( 1 ) > 0
    
    sz = size( X );
    X = [ zeros( szp( 1 ), sz( 2 ) ); X;  zeros( szp( 1 ), sz( 2 ) ) ];
    
end

% zeropad cols:
if szp( 2 ) > 0
    
    sz = size( X );
    X = [ zeros( sz( 1 ), szp( 2 )), X,  zeros( sz( 1 ), szp( 2 )) ];

end