function [u,s,v]=tensorsvd(T,ll,rl,n,absorb,cut)
%[u,s,v]=tensorsvd(T,ll,rl,n,absorb,cut)
%do svd for a tensor
%T: tensor
%ll: left legs
%rl: right legs
%n: numer of singular values: use inf for all
%absorb: 'l': left absorb
%'m': split absorb (sqrt of singular values)
%'r': right absorb
%'b': absorb singular values on both sides
%'n': no absorption of singular values
%cut: cutoff for smallest SVDs
%structure of output:
%order of legs u: leftlegs, single leg
%order of legs v: rightlegs, single leg
%CONVENTION: the conjugate of v is taken such that the original tensor is obtained
%by contracting u.s.v without taking the conmplex conjugate

[Tm,dl,dr]=lreshape(T,ll,rl);

%default: all singular values
if(nargin<4)
    n=inf;
end

%default: no absorption
if(nargin<5)
    absorb='n';
end

%default cutoff
if(nargin<6)
    cut=0;
end

% DO SVD   
try
    [u,s,v]=svd(Tm,'econ');
catch me
    disp(me.message);
    disp('ERROR OCCURED IN SVD, TRY SVDS instead:')
    [u,s,v]=svds(Tm,n);
end
       
n=min(n,size(s,1));

%check cutoff: truncate very small singular values
ds=diag(s);
ds=ds(1:n);
icut=find(ds/max(ds)<cut,1);
if(~isempty(icut))
    n=icut-1;
end

u=u(:,1:n);
v=v(:,1:n);
s=s(1:n,1:n);

% absorb singular values?
if(isequal(absorb,'m'))
   sq=diag(sqrt(diag(s)));
   u=u*sq;
   v=v*sq;
elseif(isequal(absorb,'l'))
   u=u*s;
elseif(isequal(absorb,'b'))
   u=u*s;
   v=v*s;    
elseif(isequal(absorb,'r'))    
   v=v*s;
elseif(isequal(absorb,'n'))    
   %nothing 
else
   error('undefined absorption mode');
end

ss=size(u);
u=reshape(u,[dl,ss(2)]);
v=reshape(v,[dr,ss(2)]); 

v=conj(v);

if(nargout==2)
    s=v;
end

end

