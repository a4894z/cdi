function expt = ptycho2DTPA_saveresults( sol, expt, ii )

    fprintf('\n========================================================================================================'); 
    fprintf('\nSAVING Simulated Phase Retreival Experiment, \n2D Transmission Geometry and Projection Approx Assumed...'); 

    save( expt.paths.rsdata, '*' );

    %========
    
    probe  = sol.probe; 
    sample = sol.sample; 
    spos   = sol.spos; 
    
    if isfield( sol, 'metrics' ), metrics = sol.metrics; end
    if isfield( sol, 'timings' ), timings = sol.timings; end

    A = sprintf('_%s', datestr( now, 'ddmmmyyyy_tHHMMSS' ));
    B = num2str( sol.it.epoch, '_it%d.mat');
    
    expt.paths.most_recent_date_time_save{ ii } = [ expt.paths.rsdata( 1 : end - 4 ), A, B ];
    
    try
        
        save( expt.paths.most_recent_date_time_save{ ii }, 'probe', 'sample', 'spos', 'metrics', 'timings' );
        
    catch
        
        save( expt.paths.most_recent_date_time_save{ ii }, 'probe', 'sample', 'spos' );
        
    end

    fprintf('done saving data!\n'); 
    fprintf('========================================================================================================\n\n'); 

end
