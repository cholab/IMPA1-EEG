function [F,U,i,UF]=occ(list,expected)
% Occurences (F) of unique items (U) - or desired set of 'expected' items -
% in a vector or 2D matrix 'list':
%             [F,U,i,UF]=occ(list,expected)
%and:
%  'i' is a matrix that provides indexing for each items (columns are items)
%  'UF' just puts F and U in a matrix for quick visualization
% Notes for "list'
%class support: cell, char and double
%for matrices: unique occurences are defined as unique rows


if ~nargin
    help occ
    return
end


if size(list,1)==1 && size(list,2)>1   
 list=list';
end   
no=size(list,1);


revert2double=0;
if isnumeric(list)
if isvector(list)
    [U,~,i0]=unique(list);
    if any(isnan(U))
        nan1=find(isnan(U),1);
        U=U(1:nan1);
        i0(i0>nan1)=nan1;
    end
else
   if any(isnan(list(:)))
      revert2double=1; 
      list=num2str(list);
   else
     [U,~,i0]=unique(list,'rows');
   end
end
end

switch class(list)    
    case 'char'
    [U,~,i0]=unique(list,'rows');
    case 'cell'
        try
        [U,~,i0]=unique(list);
        catch 
            k=cell(size(list));
            for n=1:numel(k)
                k{n}=class(list{n});
            end
            if all(strcmpi(k,'cell')) %nested cells
            for n=1:numel(k)
                list(n)=list{n};
                k{n}=class(list{n});
            end
            end
            if ~ all(strcmpi(k,'char')) || ~ all(strcmpi(k,'double'))
                num=find(strcmpi(k,'double'));
                 for n=1:numel(num)
                 list{num(n)}=num2str(list{num(n)});
                 end   
            end
             [U,~,i0]=unique(list);
        end
end

if exist('expected','var') && ~isempty(expected) 
  if size(expected,1)==1 && size(expected,2)>1   
    expected=expected';
  end    
  if revert2double
      expected=num2str(expected);
  end
  [~,i]=makeas(U,expected);
  [~,i0]=ismember(i0,i);
  U=expected;
end


ni=size(U,1);
i=zeros(no,ni);
for n=1:ni
   i(:,n)=n==i0; 
end

F=sum(i)';
i=logical(i);

if revert2double
    U=str2mat(U);
end

if nargout>3 
    if isnumeric(U)
      UF=[U F];
    else
      UF=[(1:numel(F))' F];
    end
    bar(F)
    set(gca,'Xtick',1:numel(F),'Xticklabel',U)
end


