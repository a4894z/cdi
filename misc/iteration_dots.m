function make_iteration_dots( alg_name, ii, ti, io )

fprintf('.');
%if (mod(ti,100)==0) || (mod(ii,100)==0), fprintf( [num2str(ti,' %d,'),' ',upper(alg_name),' ']), end
if (mod(ii,io) == 0), fprintf(num2str(ii,'<<%d>>')); end
if (mod(ti,100)==0) || (mod(ii,100)==0), fprintf( [' ',upper(alg_name),' ']), end
if mod(ti,100)==0, fprintf('\n'), end
