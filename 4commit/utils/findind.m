function [r,i]=findind(x,y,mode,constrain)
%               [r,i]=findind(x,y,'r',constrain)
%               [i,r]=findind(x,y,'i',constrain)
%It finds the closest neighbor of each x in y, and returns 
%-the corresponding ratio of the total length 'r'
%-its index 'i'
%For instance to finds the time points tp of interests in a time line tl:
%                 r=findind(tp,tl)    ['r']
%Use 'constrain' to find the closest antecedent ('>') or subsequent ('<').

if ~nargin
    help findind
    return
end

if ~exist('mode','var') || isempty(mode)
    mode='r';
end
if ~exist('constrain','var') || isempty(constrain) ||  ~(strcmp(constrain,'>') || strcmp(constrain,'<'))
   constrain=0; 
end


y=y(:);
x=x(:)';
s=length(y);
n=length(x);
x=repmat(x,[s 1]);
y=repmat(y,[1 n]);

if ~constrain
    i=abs(x-y);
elseif constrain=='>'
    i=x-y; 
    i(i<0)=NaN;
elseif constrain=='<'
    i=y-x;  
    i(i<0)=NaN;    
end
[v,i]=min(i);
i(isnan(v))=NaN;

r=i/s;
if strcmpi(mode(1),'i')
    [r,i]=swap(r,i); 
end


