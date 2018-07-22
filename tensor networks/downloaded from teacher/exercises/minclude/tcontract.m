function Q=tcontract(A,B,ma,mb,perm)
% perform a tensor contraction between indices ma of A and mb of B
% perm is an optional permutation at the end

sa=size2(A,ma);
na=numel(sa);
Ia=1:na;
Ia(ma)=[];

sb=size2(B,mb);
nb=size(sb,2);
Ib=1:nb;
Ib(mb)=[];

if sa(ma) ~= sb(mb) 
    'index dimensions do not match'
    hh
end

A=permute(A,[Ia ma]);
B=permute(B,[Ib mb]);

Q=reshape(A,[prod(sa(Ia)) prod(sa(ma))])*reshape(B,[prod(sb(Ib)) prod(sb(mb))]).';

if ~isempty(Ia) || ~isempty(Ib) 
    reexpand = [sa(Ia) sb(Ib)];    
    Q=reshape(Q,[reexpand 1]);
end

% do permutation in the end
if nargin>4
    Q=permute(Q,perm);
end

end

function mA = size2(A,ma)
% return sizes in case of D=1 legs at the end
sA = size(A);
nsA = numel(sA);
si = max(nsA,max(ma));
mA = [sA,ones(1,si-nsA)];
end