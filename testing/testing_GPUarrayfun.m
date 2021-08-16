%% Improve Performance of Element-wise MATLAB(R) Functions on the GPU using ARRAYFUN
% This example shows how |arrayfun| can be used to run a MATLAB(R) function
% natively on the GPU. When the MATLAB function contains many element-wise
% operations, |arrayfun| can provide improved performance when compared to
% simply executing the MATLAB function directly on the GPU with gpuArray
% input data. The MATLAB function can be in its own file or can be a nested
% or anonymous function. It must contain only scalar operations and
% arithmetic.

% Copyright 2010-2017 The MathWorks, Inc.

%%
% We put the example into a function to allow nested functions:
function testing_GPUarrayfun


%% Using Horner's Rule to Calculate Exponentials
% Horner's rule allows the efficient evaluation of power series expansions.
% We will use it to calculate the first 10 terms of the power series
% expansion for the exponential function |exp|. We can implement this as a
% MATLAB function.

function y = horner(x)
%HORNER - series expansion for exp(x) using Horner's rule
y = 1 + x.*(1 + x.*((1 + x.*((1 + ...
        x.*((1 + x.*((1 + x.*((1 + x.*((1 + ...
        x.*((1 + x./9)./8))./7))./6))./5))./4))./3))./2));
end

%% Preparing |horner| for the GPU
% To run this function on the GPU with minimal code changes, we could pass
% a |gpuArray| object as input to the |horner| function.  Since |horner|
% contains only individual element-wise operations, we might not realize
% very good performance on the GPU when performing each operation one at a
% time.  However, we can improve the performance by executing all of the
% element-wise operations in the |horner| function at one time using
% |arrayfun|.
%
% To run this function on the GPU using |arrayfun|, we use a handle to the
% |horner| function.  |horner| automatically adapts to different size and
% type inputs. We can compare the results computed on the GPU using both
% |gpuArray| objects and |arrayfun| with standard MATLAB CPU execution
% simply by evaluating the function directly.

hornerFcn = @horner;

%% Create the Input Data
% We create some inputs of different types and sizes, and use |gpuArray| to
% send them to the GPU.

data1  = rand( 8000, 'single' );
data2  = rand( 8000, 'double' );
gdata1 = gpuArray( data1 );
gdata2 = gpuArray( data2 );

%% Evaluate |horner| on the GPU 
% To evaluate the |horner| function on the GPU, we have two choices.  With
% minimal code changes we can evaluate the original function on the GPU by
% providing a |gpuArray| object as input. However, to improve the
% performance on the GPU call |arrayfun|, using the same calling convention
% as the original MATLAB function.
%
% We can compare the accuracy of the results by evaluating the original
% function directly in MATLAB on the CPU. We expect some slight numerical
% differences because the floating-point arithmetic on the GPU does not
% precisely match the arithmetic performed on the CPU.

gresult1 = arrayfun( hornerFcn, gdata1 );
gresult2 = arrayfun( hornerFcn, gdata2 );

comparesingle = max( max( abs( gresult1 - horner( data1 ) ) ) );
comparedouble = max( max( abs( gresult2 - horner( data2 ) ) ) );
%%
fprintf( 'Maximum discrepancy for single precision: %g\n', comparesingle );
fprintf( 'Maximum discrepancy for double precision: %g\n', comparedouble );

%% Comparing Performance between GPU and CPU
% We can compare the performance of the GPU versions to the native MATLAB
% CPU version. Current generation GPUs have much better performance in
% single precision, so we compare that.

% CPU execution
tic
hornerFcn( data1 );
tcpu = toc;

% GPU execution using only gpuArray objects 
tgpuObject = gputimeit(@() hornerFcn(gdata1));

% GPU execution using gpuArray objects with arrayfun
tgpuArrayfun = gputimeit(@() arrayfun(hornerFcn, gdata1));


fprintf( 'Speed-up achieved using gpuArray objects only: %g\n', tcpu / tgpuObject );
fprintf( 'Speed-up achieved using gpuArray objects with arrayfun: %g\n', tcpu / tgpuArrayfun );


end
