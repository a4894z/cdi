
%{

clear;close all; testing_parallel_MPI

%}

%==================================================================================================

% delete(gcp('nocreate'))
% pmode start;
% parpool(8);

% myCluster = parcluster('local');
% myCluster.NumWorkers = 10;  % 'Modified' property now TRUE
% saveProfile(myCluster);    % 'local' profile now updated,

%=======================

n = 50;
A = 2048;
float_type = 'single';
B = rand( A, A, n, float_type );
Bgpu = gpuArray( B );

%=======================
                          
a = zeros( 1, n, float_type );                   
tic
for i = 1 : n
    
    a(i) = max( abs( eig( fft2( B( :, :, i )))));
    
end
ts = toc;

%=======================

a = gpuArray( zeros( 1, n, float_type ));
tic
for i = 1 : n
    
    a(i) = max( abs( eig( fft2( gpuArray( B( :, :, i ))))));
    
end
tp = toc;

ts / tp

%=======================

a = gpuArray( zeros( 1, n, float_type ));
tic
for i = 1 : n
    
    a(i) = max( abs( eig( fft2( Bgpu( :, :, i )))));
    
end
tp = toc;

ts / tp

%=======================

a = zeros( 1, n, float_type );
tic
parfor i = 1 : n

    a(i) = max( abs( eig( fft2( B( :, :, i )))));
    
end
tp = toc;

ts / tp

%=======================

a = gpuArray( zeros( 1, n, float_type ));
tic
parfor i = 1 : n

    a(i) = max( abs( eig( fft2( Bgpu( :, :, i )))));
    
end
tp = toc;

ts / tp

%=======================

% delete(gcp('nocreate'))

%=======================

% parpool(10);
% tic
% n = 50;
% A = 512;
% a = zeros(1,n);
% parfor ( i = 1:n, 10 )
%     a(i) = max(abs(eig(rand(A,A))));
% end
% toc
% delete(gcp('nocreate'))

return

%==================================================================================================

parpool(4);

spmd
    
  if labindex == 3, labindex, end
  
  q = magic(labindex + 2);
  
end

figure
subplot( [ 1, 4, 1 ] ), imagesc(q{1});
subplot( [ 1, 4, 2 ] ), imagesc(q{2});
subplot( [ 1, 4, 3 ] ), imagesc(q{3});
subplot( [ 1, 4, 4 ] ), imagesc(q{4});


delete(gcp);

%==================================================================================================

delete(gcp('nocreate'))

parpool(4);

clear;close all;
N = 4096;

spmd
    
    codistr = codistributor1d( codistributor1d.unsetDimension, codistributor1d.unsetPartition, [ N, N ] );
                
    myLocalSize = [ N, N ]; % start with full size on each lab
    
    codistr.Dimension
    
    % then set myLocalSize to default part of whole array:
    myLocalSize( codistr.Dimension ) = codistr.Partition( labindex );
    myLocalPart = fft2(labindex * rand( myLocalSize )); % arbitrary values
    % figure( labindex ); imagesc( gather(myLocalPart) )
    D = codistributed.build( myLocalPart, codistr );
    
end

Dg = gather( D );
figure; imagesc( log10(1+abs( Dg )))
figure; spy( D == 2 );










N = 1000;
spmd
    codistr = codistributor1d(1); % 1st dimension (rows)
    C = labindex * ones(N,codistr);
end
Cg = gather( C );
figure; imagesc( Cg )






clear;
N = 5;
Xs = magic(N);          % Replicated on every worker
spmd

    X = rand(N);          % Replicated on every worker
    C1 = codistributed( X ); % Partitioned among the workers
    
    if labindex == 3

        X
        
    end
    
    if labindex == 2

        X
    end
    
end

Xs2 = X{2};








spmd
    
    N = 10;
    X = rand(N);          % Replicated on every worker
    C1 = codistributed(X); % Partitioned among the workers
    if labindex == 1
       C1(LocalPart )
    end
end

if labindex == 1
    Z = gather( C1 )
    X
end





clear; close all;

delete(gcp('nocreate'))

parpool(4);

spmd
    N = 1001;
    globalSize = [N,N];
    % Distribute the matrix over the second dimension (columns),
    % and let the codistributor derive the partition from the 
    % global size.
    codistr = codistributor1d(2, ...W
                 codistributor1d.unsetPartition,globalSize)
 
    % On 4 workers, codistr.Partition equals [251,250,250,250].
    % Allocate storage for the local part.
    localSize = [ N, codistr.Partition( labindex ) ]
    L = zeros(localSize);
    
    % Use globalIndices to map the indices of the columns 
    % of the local part into the global column indices.
    % On 4 workers, globalInd has the values:
    % 1:251    on worker 1
    % 252:501  on worker 2
    % 502:751  on worker 3
    % 752:1001 on worker 4
    globalInd = codistr.globalIndices( 2 ); 

    
    % Initialize the columns of the local part to the correct value.
    for localCol = 1 : length( globalInd )
        
        globalCol = globalInd( localCol );
        L( :, localCol ) = globalCol;
        
    end
    
    D = codistributed.build( L, codistr );
end

figure; imagesc(D)


x1 = [localSize{:}];
x2 = [globalInd{1}];


[codistr{1}]


help codistributor1d













spmd
%     C = zeros( 2, 22,codistributor1d( 2, [ 6, 6, 5, 5 ] ));
    C = zeros( 2, 21, codistributor1d( 2 ));
    
    if labindex == 1
       K = globalIndices( C, 2 )     % returns K = 1:6.
       [ E, F ] = globalIndices( C, 2 )
    end
    
    if labindex == 2
        K = globalIndices( C, 2 )
       [ E, F ] = globalIndices( C, 2 ) % returns E = 7, F = 12.
    end
    if labindex == 3
        K = globalIndices( C, 2 )
       [ E, F ] = globalIndices( C, 2 ) % returns E = 7, F = 12.
    end
    if labindex == 4
        K = globalIndices( C, 2 )
       [ E, F ] = globalIndices( C, 2 ) % returns E = 7, F = 12.
    end
    K = globalIndices( C, 2, 3 );      % returns K = 13:17.
    [ E, F ] = globalIndices( C, 2, 4 );  % returns E = 18, F = 22.
    
 end











N = 512;
A = rand( N, N );

parfor k = 1 : 200
   a( k ) = fft( A( :, kk ), [], 2 );
end

for k = 1 : 200
   a( k ) = fft( A( :, kk ), [], 2 );
end

parfor k = 1 : 200
    if labindex == 3, labindex, end
end
 

B = zeros( N, N );


spmd
  LP = getLocalPart(A);
  LPf = fft(LP);
  T = codistributed.build(LPf, getCodistributor(A));
end





spmd 
    
    if labindex == 3
        labindex, 
    end
    
    a = fft( A( :, labindex ), [], 2 );
    
    if labindex == 3
        labindex, 
        a
    end
    
    B( :, labindex ) = a;
    
end


































delete(gcp('nocreate'))



%==================================================================================================



















 M = 0;                     % M specifies maximum number of workers
 y = ones(1,100);
 parfor ( i = 1:100, M )
      y(i) = i;
 end
 
 

cluster = parcluster;

values = [3 3 3 7 3 3 3];
parfor (i=1:numel(values),cluster)
    out(i) = norm(pinv(rand(values(i)*1e3)));
end

 
 
 
 
 
delete(gcp('nocreate'))

tic
ticBytes(gcp);
n = 200;
A = 500;
a = zeros(1,n);
parfor i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
tocBytes(gcp)
toc