function [X,other]=ref(X,reference, forcelast, existingref)
%        [X,other]=ref(X,reference, forcelast, existingref)
%It checks reference status, add the reference if missing and rereference
%data to "reference"  if provided (reference=0 for average rereferencing).
%Any 'other' flat channel is neglected and returned.


if ~nargin
    help ref
    return
end
if ~exist('forcelast','var') || isempty(forcelast)
    forcelast=1;
end


el_in=size(X,1);  other=0;
%% CHECK REFERENCING STATUS    
ref_in=find(all(X(:,:,1)==0,2));
if numel(ref_in)>1
    other=ref_in(1:end-1);
    ref_in=ref_in(end);
    display(['Neglecting flat channels ' num2str(other(:)')])      
end
if forcelast && ~isempty(ref_in) && ref_in~=el_in
   display(['Forcing to neglect flat channel ' num2str(ref_in)])  
   other=ref_in;
   ref_in=[]; 
end
if ~isempty(ref_in)
    display(['Assuming reference is channel ' num2str(ref_in)])        
elseif max(abs(nanmean(X(:,:))))<mean(max(abs((X(:,:)))))/10000
    ref_in=0;
    display('Data appears to be average rereferenced.')        
else   
    el_in=el_in+1;
    ref_in=el_in;
    display(['Assuming reference is missing channel ' num2str(ref_in) ' (putting it back)'])        
    X(ref_in,:,:)=0; 
%     warning('Reference status could not be determined!')                   
%     display('This could mean channels were already removed, "applymask" exits')                           
end
if ~exist('reference','var') || isempty(reference)
    reference=ref_in;
end

if ref_in~=reference
    if ~reference %average rereference
       fprintf('Average re-referencing...') 
       if any(isnan(X(:)))
       X=X-meandim(X,1,[],2);
       else
       M=eye(el_in)-ones(el_in)*1/el_in;
       X=reshape(M*X(:,:),size(X,1),size(X,2),[]);
       end
       fprintf(' done.\n')        
    else
       X=reref(X,reference);
    end 
end


