function [ TFview ] = getview_2DsampleTF( TF, vs_r, vs_c, spos, shifttype )

if ~exist( 'shifttype', 'var' )
  
    shifttype = 'px'; 
    
end

%=======================

if strcmp( shifttype, 'subpx' )
    
    TFb = subpixelshift2D( TF, spos );
%     TFview = TFb( vs( 1, : ), vs( 2, : ));     
    TFview = TFb( vs_r, vs_c ); 
    
else

%     TFb = nocircshift2D( TF, round( spos ));
%     TFview = TFb( vs.r, vs.c );
    
%     TFb = circshift( TF, round( spos ));
%     TFview = TFb( vs.r, vs.c );

%     TFview = TF( round( vs( 1, : ) - 1.0 * spos( 1 )), round( vs( 2, : ) - 1.0 * spos( 2 )));
    TFview = TF( round( vs_r - 1.0 * spos( 1 )), round( vs_c - 1.0 * spos( 2 )));
    
end   

%=======================  








%{


A = [ 5, 4, 2, 6, 5; ...
      1, 5, 6, 2, 8;
      17, 4, 9, 1, 11 ]; 
  
A = uint32( A );
  
B = (1 : 20)'

C = B( A )































A = [0.1299    0.3371    0.5285; ...
     0.5688    0.1622    0.1656; ...
     0.4694    0.7943    0.6020; ...
     0.0119    0.3112    0.2630];
 
B = [1 2 4];
C = [2 3 4];

D = zeros(size(A));
for k = 1:numel(B)
  D(B(k):C(k), k) = A(B(k):C(k), k);
end






[s1, s2] = size(A);
sp = [s1 + 1, s2];
M = zeros(sp);
M(sub2ind(sp, B, 1:s2)) = 1;
M(sub2ind(sp, C + 1, 1:s2)) = -1;
M = cumsum(M(1:s1, :), 1);
D = A .* M;





A  = rand(1000, 1000);
BB = randi(1000, 1, 1000);
CC = randi(1000, 1, 1000);
B  = min(BB, CC);
C  = max(BB, CC);


tic; 

for q = 1:100
    
D = zeros(size(A));
for k = 1:numel(B)
  D(B(k):C(k), k) = A(B(k):C(k), k);
end

end

toc






tic; 

for q = 1:100
    

    [s1, s2] = size(A);
    sp = [s1 + 1, s2];
    M = zeros(sp);
    M(sub2ind(sp, B, 1:s2)) = 1;
    M(sub2ind(sp, C + 1, 1:s2)) = -1;
    M = cumsum(M(1:s1, :), 1);
    D = A .* M;

end

toc







%}

