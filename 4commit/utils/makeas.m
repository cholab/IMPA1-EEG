function [X,i]=makeas(keyA,keyB,A,filler)
%        [keyA_as_keyB,i]=makeas(keyA,keyB)
%        [A_as_keyB,        i]=makeas(keyA,keyB,A)
%It makes keyA as keyB, i.e. same size and same order of matching items.
%'i' can be used to reorganize any list that matches keyA dimensions, but
%note that missing items in keyA are NaNs in 'i'
%Optionally, the ordering scheme can be applied directly to A.
%Keys are numeric/char or cells of strings (numeric is faster), A can be anything.


if ~nargin
    help makeas
    return
end


if exist('A','var') && ~isempty(A) 
    if length(A)== length(keyB)
        display('**  Input/key missmatch: keys order might be reversed!! **')        
         [keyA,keyB]=swap(keyA,keyB);
    elseif length(A)~= length(keyA)
        error('Unsolvable input/key missmatch!')
    end
else
    A=keyA;
end

[keyA,keyB]=recode(keyA,keyB);
[present,i]=ismember(keyB,keyA); 
i(~present)=NaN;
X(1:numel(i),1:size(A,2))=A(1,1); 
X(present,:)=A(i(present),:);    
if ~all(present)
%use some fillers for missing items  
if exist('filler','var')
    if strcmp(class(filler),class(A(1)))
    if ismember(filler,A)
        error('Filler matches input!')
    end
    else
        warning('Filler-input class missmatch! - default filler is used')
        clear filler
    end
end

if ~exist('filler','var')
if ischar(A)
    filler=' ';
elseif iscell(A)
    filler=cell(1);
else
    filler=NaN;
end
end

  X(~present,:)=filler;  
end
 

