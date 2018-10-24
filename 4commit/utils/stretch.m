function [X]=stretch(X,permute_order);
% Stretch converts the multidimensional matrix X into a 2d matrix. Each
% higher dimension of X is vectorized, maintaining row only. The
% permute_oder argument allows you to reshape X to choose the organizing
% row.


if exist('permute_order','var')&&~isempty(permute_order)
    X=permute(X,permute_order);
end

X=X(:,:);
    
end


