function [ sol, expt ] = define_initial_misc( sol, expt )

    sol.csys   = expt.csys;
    sol.lambda = expt.lambda;
    sol.sz     = expt.sz;
    
    %=============================================
    % Performance metric and epoch update counters
    %=============================================
    
    if ~isfield( sol, 'it' ),       sol.it = struct;  end
    
    if ~isfield( sol.it, 'mtot' ),  sol.it.mtot  = 1; end
    if ~isfield( sol.it, 'metr' ),  sol.it.metr  = 1; end
    if ~isfield( sol.it, 'epoch' ), sol.it.epoch = 1; end
    
    %=====================================================
    % Modify blemish array that defines unmeasured regions
    %=====================================================
    
    
%     tmp0 = fftshift( expt.meas.blemish );
%     
% %     figure; imagesc(tmp0)
% 
%     tmp0( :, 561 : 563 ) = 1;
%     
%     expt.meas.blemish = fftshift( tmp0 );
    
    
    %==========================================
    % fliplr/flipud/rot90()/etc of measurements
    %==========================================

%     expt.meas.D = flip( expt.meas.D, 2 );
%     expt.meas.D = flip( expt.meas.D, 1 );
%     expt.meas.D = nocircshift3D( expt.meas.D, [ 0, +2, 0 ] );

%     tmp0 = sol.spos.rs;
%     sol.spos.rs( :, 1 ) = tmp0( :, 2 );
%     sol.spos.rs( :, 2 ) = tmp0( :, 1 );
%     
%     clear( 'tmp0' )

%     expt.meas.D = nocircshift3D( fftshift( fftshift( expt.meas.D, 1 ), 2 ), [ +4, 0, 0 ] );
%     expt.meas.D = fftshift( fftshift( expt.meas.D, 1 ), 2 );

end

