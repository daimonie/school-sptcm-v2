function disp0(varargin)
% simpler output of strings mixed with numbers
% example: disp0('the result is:',2, ' and ',3);

s='';   
ls=size(varargin,2);

for k=1:ls 
    if(isstr(varargin{k}));
        s=[s varargin{k}];
    else
        s=[s num2str(varargin{k})];
    end

end
disp(s);

end