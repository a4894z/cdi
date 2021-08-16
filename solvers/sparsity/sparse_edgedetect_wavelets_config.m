function [ sparsity ] = sparse_edgedetect_wavelets_config( sparsity, N )


% 'haar', 'db4', 'coif4', 'coif2'
% 'bior1.1', 'bior1.3', 'bior1.5', 'bior2.2', 'bior2.4', 'bior2.6', 'bior2.8', 'bior3.1', 'bior3.3', 'bior3.5', 'bior3.7''bior3.9', 'bior4.4', 'bior5.5', 'bior6.8'

% 'coif1', ... , 'coif5'
% 'sym2', ... , 'sym8', ... ,'sym45'
% 'fk4', 'fk6', 'fk8', 'fk14', 'fk22'


sparsity.num_levels = 3;

sparsity.type = 'haar';
%sparsity.type = 'db2';
%sparsity.type = 'db6';
%sparsity.type = 'coif2';
%sparsity.type = 'fk4';
%sparsity.type = 'bior3.5';
%sparsity.type = 'sym8';

sparsity.firm_mult = 0.5;

%=======================

% LEVEL 1 (FINEST DETAILS)

% amplitude
sparsity.level(1).A = round( 1 * N.NrNc );      
%sparsity.shrink{1}.A = 'soft';
%sparsity.shrink{1}.A = 'firm';   
sparsity.shrink{1}.A = 'hard';

% horizontal detail
sparsity.level(1).H = round( 0.2 * N.NrNc );    
sparsity.shrink{1}.H = 'soft';
%sparsity.shrink{1}.H = 'firm'; 
%sparsity.shrink{1}.H = 'hard';

% vertical detail
sparsity.level(1).V = round( 0.2 * N.NrNc );
sparsity.shrink{1}.V = 'soft';
%sparsity.shrink{1}.V = 'firm';
%sparsity.shrink{1}.V = 'hard';

% diagonal detail
sparsity.level(1).D = round( 0.2 * N.NrNc );
sparsity.shrink{1}.D = 'soft';
%sparsity.shrink{1}.D = 'firm';
%sparsity.shrink{1}.D = 'hard';

%=======================

% LEVEL 2

% amplitude
sparsity.level(2).A = round( 1 * N.NrNc );
sparsity.shrink{2}.A = 'soft';
%sparsity.shrink{2}.A = 'firm';
%sparsity.shrink{2}.A = 'hard';

% horizontal detail
sparsity.level(2).H = round( 0.2 * N.NrNc );
%sparsity.shrink{2}.H = 'soft';
%sparsity.shrink{2}.H = 'firm';
sparsity.shrink{2}.H = 'hard';

% vertical detail
sparsity.level(2).V = round( 0.2 * N.NrNc );
%sparsity.shrink{2}.V = 'soft';
%sparsity.shrink{2}.V = 'firm';
sparsity.shrink{2}.V = 'hard';

% diagonal detail
sparsity.level(2).D = round( 0.2 * N.NrNc );
%sparsity.shrink{2}.D = 'soft';
%sparsity.shrink{2}.D = 'firm';
sparsity.shrink{2}.D = 'hard';

%=======================

% LEVEL 3

% amplitude
sparsity.level(3).A = round( 1 * N.NrNc );
%sparsity.shrink{3}.A = 'soft';
%sparsity.shrink{3}.A = 'firm';
sparsity.shrink{3}.A = 'hard';

% horizontal detail
sparsity.level(3).H = round( 0.3 * N.NrNc );
%sparsity.shrink{3}.H = 'soft';
%sparsity.shrink{3}.H = 'firm';
 sparsity.shrink{3}.H = 'hard';

% vertical detail
sparsity.level(3).V = round( 0.3 * N.NrNc );
%sparsity.shrink{3}.V = 'soft';
%sparsity.shrink{3}.V = 'firm';
sparsity.shrink{3}.V = 'hard';

% diagonal detail
sparsity.level(3).D = round( 0.3 * N.NrNc );
%sparsity.shrink{3}.D = 'soft';
%sparsity.shrink{3}.D = 'firm';
sparsity.shrink{3}.D = 'hard';


%=======================

% LEVEL 4

% amplitude
sparsity.level(4).A = round( 0.3 * N.NrNc );
% sparsity.shrink{4}.A = 'soft';
sparsity.shrink{4}.A = 'firm';
% sparsity.shrink{4}.A = 'hard';

% horizontal detail
sparsity.level(4).H = round( 0.2 * N.NrNc );
% sparsity.shrink{4}.H = 'soft';
sparsity.shrink{4}.H = 'firm';
% sparsity.shrink{4}.H = 'hard';

% vertical detail
sparsity.level(4).V = round( 0.2 * N.NrNc );
% sparsity.shrink{4}.V = 'soft';
sparsity.shrink{4}.V = 'firm';
% sparsity.shrink{4}.V = 'hard';

% diagonal detail
sparsity.level(4).D = round( 0.2 * N.NrNc );
% sparsity.shrink{4}.D = 'soft';
sparsity.shrink{4}.D = 'firm';
% sparsity.shrink{4}.D = 'hard';


%=======================

% LEVEL 5

% amplitude
sparsity.level(5).A = round( 0.3 * N.NrNc );
% sparsity.shrink{5}.A = 'soft';
sparsity.shrink{5}.A = 'firm';
% sparsity.shrink{5}.A = 'hard';

% horizontal detail
sparsity.level(5).H = round( 0.2 * N.NrNc );
% sparsity.shrink{5}.H = 'soft';
sparsity.shrink{5}.H = 'firm';
% sparsity.shrink{5}.H = 'hard';

% vertical detail
sparsity.level(5).V = round( 0.2 * N.NrNc );
% sparsity.shrink{5}.V = 'soft';
sparsity.shrink{5}.V = 'firm';
% sparsity.shrink{5}.V = 'hard';

% diagonal detail
sparsity.level(5).D = round( 0.20 * N.NrNc );
% sparsity.shrink{5}.D = 'soft';
sparsity.shrink{5}.D = 'firm';
% sparsity.shrink{5}.D = 'hard';

%==================================================================================================



