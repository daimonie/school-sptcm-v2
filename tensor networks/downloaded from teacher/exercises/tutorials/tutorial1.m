%% Tutorial 1: First MATLAB warmup
% Let's try out a few things:
%%


%% Generate a random vector of length 4
% type 'help rand' to see the help function of rand
% initialize random number generator
rand('seed',1);
v = rand(4,1);
% output it: leave away the semicolon at the end of the line
v
% or use
disp('output v:');
disp(v);


%% Generate matrices 
Mid = eye(2)    % 2x2 Identity matrix
Mdiag = diag([1 2 3]) % diagonal matrix with 1 2 3 on the diagonal
M = rand(3,3)  % a random 3x3 matrix


%% Accessing parts of a matrix
M(:,1)  % print first column
M(2,:)  % print second row
M(1:2,1:2) % plot the 2x2 submatrix

%% Useful matrix operations
Mconj = M';  % matrix conjugation
a = M + 2;   % add 2 to all elements
b = M * M;   % matrix-matrix multiplication
c = M .* M;  % element-wise multiplication
d = expm(M); % matrix exponential

%% Useful matrix routines
% use 'help _functionname_' to see the syntax of these functions
M = rand(10,10);    % random 10x10 matrix
M = M + M';         % initialize hermitian matrix
[V,D] = eig(M);     % Eigenvalue decomposition
% Iterative eigensolver: compute the 2 eigenvector with smallest
% eigenvalues
[V,D] = eigs(M,2,'SA');
% Singular value decomposition (economy size)
[U,s,V] = svd(M,'econ');
% qr decomposition
[Q,R] = qr(M,0);

%% Generating tensors
T1 = rand(2,3,3);    % 2x3x3 random tensor
T2 = zeros(2,4,1);   % 2x4x1 tensor with zeros

%% Get size of a tensor
size(T1)
disp(['Total number of elements:',num2str(numel(T1))]);


%% reshaping 
v = 1:9    % this is a vector with elements from 1...9
M = reshape(v,3,3) % reshape it into a 3x3 matrix
v2 = 1:8   % this is a vector with elements from 1...8
T = reshape(v2,2,2,2) % reshape it into a 2x2x2 tensor

%% permutation
Mp = permute(M,[2 1]) % permute indices of matrix (=transpose)
Tp = permute(T,[1 3 2]) % permute indices of tensor

%% Tensor-tensor contraction using matrix-matrix multiplication
% NOTE: we will use special functions tcontract and ncon which will do
% this automatically (see tutorial 5)!
A = rand(2,3,4); 
B = rand(3,4,2);
% contract legs 2 and 3 of A with legs 1 and 2 of B
% determine sizes of A and B
sA = size(A);
lAdimleft = sA(1); 
lAdimright = sA(2:3); 
sB = size(B);
lBdimleft = sB(1:2);
lBdimright = sB(3);
% now reshape A and B into matrices
MA = reshape(A,lAdimleft,prod(lAdimright));
MB = reshape(B,prod(lBdimleft),lBdimright);
% now multiply matrices
MC = MA*MB;
% and now reshape back intro tensor
disp('Resulting 2x2 tensor:');
C = reshape(MC,[lAdimleft,lBdimright])

%% Tensor-tensor contraction requiring a permutation first
% Same example as before but now the legs are not in correct order
A = rand(3,2,4);
B = rand(2,4,3);
% contract legs 1 and 3 of A with legs 3 and 2 of B
% first permute the legs: 
% the legs of A which are going to be contracted need to be on the right
% side:
A = permute(A,[2,1,3]); % now is a 2x3x4 tensor 
B = permute(B,[3,2,1]); % now is a 3x4x2 tensor 
% and now use same code as above:
% determine sizes of A and B
sA = size(A);
lAdimleft = sA(1); 
lAdimright = sA(2:3); 
sB = size(B);
lBdimleft = sB(1:2);
lBdimright = sB(3);
% now reshape A and B into matrices
MA = reshape(A,lAdimleft,prod(lAdimright));
MB = reshape(B,prod(lBdimleft),lBdimright);
% now multiply matrices
MC = MA*MB;
% and now reshape back intro tensor
disp('Resulting 2x2 tensor:');
C = reshape(MC,[lAdimleft,lBdimright])



%% Matlab cell arrays
% Cell arrays are useful containers in Matlab, which can be filled with
% arbitrary types (also mixing types in the same cell), see e.g.:
c1 = cell(1,2); % initialize a 1x2 cell with empty elements
c1{1} = rand(2,2);  % add a random 2x2 matrix to c1
c1{2} = rand(3,3);  % add a random 3x3 matrix to c2
c1{end+1} = 2;      % add a number to the cell array
c1      % output content of cell array
c1{1}   % output the first element
c1{2}   % output the second element




