%% Tutorial 5: Additional tensor network and other tools
% Check out these useful routines, especially ncon and tensorsvd simplify
% tensor network codes a lot!
%%

%% addpath
% First you need to add the path to the additional functions tested here:
addpath('../minclude');

%% disp0
% this function just simplifies output of text/numbers:
% instead of writing
disp(strcat('The result is:',num2str(4)));
% simply write
disp0('The result is:',4);


%% lreshape
% with lreshape you can reshape and permute some legs and create a matrix
% directly. Checkout the documentation and/or code by 'help lreshape'
% Let's take the following example:
rand('seed',1);
A = rand(3,2,5,4);
B = rand(2,5,3);
% contract leg 3 of A with leg 2 of B 
% This becomes:
[MA,lAdimleft,lAdimright] = lreshape(A,[1 2 4],[3]);
[MB,lBdimleft,lBdimright] = lreshape(B,[2],[1 3]);
MC = MA*MB;
C = reshape(MC,[lAdimleft,lBdimright]);
disp0('Resulting size of C: ',size(C));


%% tcontract
% tcontract allows you to contract two tensors (same as above) in one command
C = tcontract(A,B,3,2);
disp0('Resulting size of C: ',size(C));
% another example:
% contract A and B on indices [1,2] of A and [3,1] of B
D = tcontract(A,B,[1 2],[3 1]);
disp0('Resulting size of D: ',size(D));


%% ncon
% ncon is a great tool to contract entire tensor networks, written by
% Robert N. C. Pfeifer, Glen Evenbly, Sukhwinder Singh, Guifre Vidal
% see documentation https://arxiv.org/abs/1402.0939
% With this function above example is done in one line:
A = rand(3,2,5,4);
B = rand(2,5,3);
C = ncon({A,B},{[-1 -2 1 -3],[-4 1 -5]});
disp0('Resulting size of C: ',size(C));


%% Example: Doing an SVD of a tensor with lreshape
T = rand(1,2,3,2);   
% lets to an svd of this tensor with legs 1 and 3 on the left
% and legs 2 and 4 on the right:
[MT,dl,dr] = lreshape(T,[1 3],[2 4]);
[MU,s,MV] = svd(MT,'econ');
disp('Singular values:');
disp(diag(s)')
% now reshape U and V back to tensors:
U = reshape(MU,[dl,size(s,1)]);
V = reshape(MV,[dr,size(s,2)]);
disp0('Size of U: ',size(U));
disp0('Size of S: ',size(s));
disp0('Size of V: ',size(V));

%% tensorsvd
% this function does above in one line:
[U,s,V] = tensorsvd(T,[1 3],[2 4]);
disp0('Size of U: ',size(U));
disp0('Size of S: ',size(s));
disp0('Size of V: ',size(V));
% WATCH OUT: in the usual svd routine the V contains the complex conjugates
% whereas tensorsvd does NOT return the complex conjugates
% In other words: By contracting U.s.V one recovers the matrix T


%% tensorsvd with compression
% One can specify the number of singular values to keep:
nsingularvalues=3;
[U,s,V] = tensorsvd(T,[1 3],[2 4],nsingularvalues);
disp0('Size of U: ',size(U));
disp0('Size of S: ',size(s));
disp0('Size of V: ',size(V));


%% Doing an SVD of a tensor and multipy resulting tensors back together
% Check out the following example and write down the diagrams by hand
T = rand(1,2,3,2);  
[U,s,V] = tensorsvd(T,[1 3],[2 4]);
T2 = ncon({U,s,V},{[-1 -3 1],[1 2],[-2 -4 2]});
% compute the sum of the abs value of all elements:
sumabsdiff = sum(abs(T(1:end)-T2(1:end)));
disp0('Difference between T and T2:',sumabsdiff);


%% Absorbing singular values in tensorsvd
% With tensorsvd you can optionally absorb the singular values in U
% (with 'l') or in V ('r') or sqrt(s) on both sides ('m'):
T = rand(1,2,3,2);  
[Us,s,V] = tensorsvd(T,[1 3],[2 4],inf,'l'); %absorb left
%singular values are contained in Us 
T2 = ncon({Us,V},{[-1 -3 1],[-2 -4 1]});
% compute the sum of the abs value of all elements:
sumabsdiff = sum(abs(T(1:end)-T2(1:end)));
disp0('Difference between T and T2:',sumabsdiff);


