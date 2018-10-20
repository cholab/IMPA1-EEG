function [out,out2]=convertind(in,dims)
%Indicing strategy converter: it inteprets/check input and changes
%subscripts into logical format and vicersa.
%If 'dims' is NOT provided:
%                     sub=convertind(logical);    logical to linear indexing
%           [row,column]=convertind(logical);      logical to subscripts
%If 'dims' is provided:
%                   lind=convertind(sub,dims);
%                   lind=convertind([row column],dims);
% 
%Note that:
%-matrix dimension dims is needed for subscript to logic indexing
%-for matrices subscripts are different from linear indexing, here
%subscripts are returned if two output arguments are required i.e [rows columns]
%-numerical input is supported and handled as *binary*, logical indexing (i.e. in=in~=0);

if~nargin
    help convertind
    return
end

if exist('dims','var')
    if numel(dims)==1
        dims=[dims 1];
    end
    if ~isvector(in)
       in=sub2ind(dims,in(:,1),in(:,2));
    end
    out=false(dims);
    out(in)=1;
    
else %logic or numerical binary form is assumed
    if ~islogical(in)
       in=logical(in);
    end
    if nargout==2
      [out,out2]=find(in);
    else
      out=find(in);
    end
end
