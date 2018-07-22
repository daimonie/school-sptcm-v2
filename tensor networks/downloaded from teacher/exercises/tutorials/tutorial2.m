%% Tutorial 2: Note on functions in Matlab
% Functions are usually definied in separate .m files
% It is not possible to combine function definitions within scripts.
% It is possible to add several functions within the same file as below,
% but then only the top function is accessible from external scripts

%% Define a function
% The name of the m-file needs to be the same as the function, and then you 
% can call it from other scripts or from the command line.
% E.g. type 'tutorial(1,3)' on the command line to call this function
% You can also define functions which do not return anything
function res = tutorial2(a, b)
    % check number of input arguments with keyword nargin
    if nargin == 0
        % no input: specify default values:
        a = 1;
        b = 1;
    end
    % Call another function defined below (not accessible from outside)
    res = myfunc(a) + myfunc(b);
end

function y = myfunc(x)
    % Example function which can be called from tutorial2 function
    y = x^3 + x;
end