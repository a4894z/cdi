function [ scpm ] = compute_scpm_photonocc( probe )

scpm.fro2    = squeeze( sum( sum( abs( probe ) .^ 2 )));
scpm.fro2TOT = sum( scpm.fro2 );
scpm.occ     = scpm.fro2 / scpm.fro2TOT;            % INTENSITY IS A MEASURE OF NUMBER OF PHOTONS, NOT AMPLITUDE 
scpm.N       = length( scpm.occ );

