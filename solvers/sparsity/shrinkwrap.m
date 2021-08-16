function [ V, Rthresh ] = shrinkwrap( V, swparams )    

%{

if ~isfield( swparams, 'blurx' ), swparams.blurx = 0.02; end
if ~isfield( swparams, 'blury' ), swparams.blury = 0.02; end

if ~isfield( swparams, 'sparselvl' )
%     warning('MUST SPECIFY SHRINKWRAP SPARSITY LEVEL, EXITING.'); 
    warning('MUST SPECIFY SHRINKWRAP SPARSITY LEVEL, SETTING TO 50%');  
    swparams.sparselvl = 0.50;
end

%}

sz = single( size( V ));
% sz = uint64( sz );

% gaussian blur probe, disallow negative values:
W = abs( V );
lpfBLUR = single( [ sz( 1 ) * swparams.blury, sz( 2 ) * swparams.blurx ] );
[ Wb ] = lpf_gauss( W, lpfBLUR );
Wb( Wb < 0 ) = 0;

tau = find_thresh_from_sparsitylevel( Wb, round( sz( 1 ) * sz( 2 ) * swparams.sparselvl ) );
% tau = swparams.pct_of_max * max( Wb( : ));

%     [ V, Rthresh ] = soft_shrinkage( V, W, Wb, tau );
%     [ V, Rthresh ] = firm_shrinkage( V, W, Wb, 0.5 * tau, tau );
[ V, Rthresh ] = hard_shrinkage( V, Wb, tau );