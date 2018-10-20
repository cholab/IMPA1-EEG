function Y=averager(X,dim,pcabased)
%                         Y=averager(X,dim,pcabased)
% It averages 'X' (or the reconostructed signal after retaing the first component if 'pcabased') along dimensions 'dim'.
% If 'pcabased' is a number or a vector, instead of the averages the value for the corresponding sources is kept.
%NOTE: 3D support for pcabased will extract components iteratively through
%the 1d e.g. if rows are frequencies and columns are timepoints, components
%will explain variability for each frequency indipendently, this is useful
%to address var(1/f) relations

if ~nargin
    help averager
else
if ~exist('pcabased','var') || isempty(pcabased) || all(pcabased==0)
    Y=meandim(X,dim);
else
    
    if islogical(pcabased) && pcabased==1
        avmode=true;
    else
        avmode=false;
    end
    d=size(X);   
    if numel(d)==2        
       Y=cfilter(X);
    if avmode
    Y=meandim(Y,dim);
    else
    if  dim==1       
    Y=Y(pcabased,:);        
    else
    Y=Y(:,pcabased);                
    end
    end
               
    else
    Y=zeros(d(1:2));
    for n=1:d(1) 
        y=cfilter(squeeze(X(n,:,:)));
    if avmode
        Y(n,:)=meandim(y,2);
    else
        Y(n,:)=y(:,pcabased);
    end
    end
    end
    

end
end