%

%{

cd /net/s8iddata/export/8-id-ECA/Analysis/atripath/cdi

clear; close all; runsolver_ptycho2DTPA_multiple_trials

%}

%====================================================================================================================================================

paths_run = {};

paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu1/cdi_ePIE_full_randT/';
paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu1/cdi_rPIE_full_alpha01_randT/';
paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu1/cdi_rPIE_full_alpha05_randT/';
paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu1/cdi_rPIE_full_alpha10_randT/';
paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu1/cdi_rPIE_full_alpha25_randT/';
paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu1/cdi_rPIE_full_alpha50_randT/';

%{

% test paths
for pp = 1 : length( paths_run ), 
    cd( paths_run{ pp }); 
end

%}

%========

N_trials = 10;

[ ~, ~ ] = run_trials( paths_run, N_trials );

%====================================================================================================================================================

% paths_run = {};
% 
% paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu4/cdi_rPIE_stochGD_alpha06/';
% paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu4/cdi_rPIE_stochGD_alpha07/';
% paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu4/cdi_rPIE_stochGD_alpha08/';
% paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu4/cdi_rPIE_stochGD_alpha09/';
% paths_run{ end + 1 } = '/net/s8iddata/export/8-id-ECA/Analysis/atripath/gpu4/cdi_rPIE_stochGD_alpha10/';

% % test paths
% for pp = 1 : length( paths_run ) 
%     cd( paths_run{ pp }); 
% end

% %========
% 
% N_trials = 10;
% 
% [ ~, ~ ] = run_trials( paths_run, N_trials );

%====================================================================================================================================================

function [ sol, expt ] = run_trials( paths_run, N_trials )

    for pp = 1 : length( paths_run )

        %========

        cd( paths_run{ pp })

        %========

        A = [ 'independenttrials', sprintf('_%s', datestr( now, 'ddmmmyyyy_tHHMMSS' )) ];

        mkdir( A ); 

        %========

        for tt = 1 : N_trials

            [ sol, expt ] = runsolver_ptycho2DTPA;

            %========

            B = [ A, '/', num2str( tt, 'trial_%d') ];
            mkdir( B )

            try movefile( '*.jpg', B ); catch, end

            try 

                for ii = 1 : length( expt.paths.most_recent_date_time_save )

                    movefile( expt.paths.most_recent_date_time_save{ ii }, B );

                end

            catch

            end

            try movefile( './sim_ptycho2DTPA.mat', B ); catch, end

            try copyfile( './sim_ptycho2DTPA_0.mat', './sim_ptycho2DTPA.mat' ); catch, end

        end

    end

end