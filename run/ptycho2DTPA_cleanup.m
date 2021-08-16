function [ sol, expt ] = ptycho2DTPA_cleanup( sol, expt )

    clearvars -except expt sol

    sol.phi    = []; sol = rmfield( sol, 'phi' );
    sol.phiOLD = []; sol = rmfield( sol, 'phiOLD' );
    sol.unmeas = []; sol = rmfield( sol, 'unmeas' );

end