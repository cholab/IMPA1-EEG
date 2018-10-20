function i=extendtwin(i,lr,method)
%                         i=extendtwin(i,lr,method)
%It extends time windows marked as i==1 of  'lr(1)' samples  to the left and
%'lr(2)' samples to the right, or symmetrically of 'lr' timepoints,
%according to 'method'  [binary] :
%        -'binary'        binary filling, i.e. 1 over the provided sample length
%        -'sinus'    sinusoidal filling, i.e. continous decay to 0 over the provided sample length


if ~nargin
   help extendtwin
   return
end

if numel(lr)==1
    lr=[lr lr];
end
if ~exist('method','var') || isempty(method)
    method='binary';
end

if size(i,2)<size(i,1)
    i=i';
end

i=[zeros(1,lr(1)) 0 i==1 0 zeros(1,lr(2))];
iA=[i(2:end) 0];
iZ=[0 i(1:end-1)];

A=find(~i & iA)+1;
Z=find(~i & iZ)-1;

A=repmat((-lr(1):-1), [numel(A) 1] )+ repmat(A', [1 lr(1)] );
Z=repmat((1:lr(2)),   [numel(Z) 1] )+   repmat(Z', [1 lr(2)] );

if strcmpi(method(1),'s')
    s=0.5+sin(linspace(-pi/2,pi/2,lr(1)+2))/2;
    iA=i;
    iA(A)=repmat(s(2:end-1),[size(A,1) 1]);
    
    s=fliplr(0.5+sin(linspace(-pi/2,pi/2,lr(2)+2))/2);    
    iZ=i;
    iZ(Z)=repmat(s(2:end-1),[size(Z,1) 1]);
    i= max([i ;iA; iZ]);
else
    A=[A(:) ;Z(:)]; 
    A(A<1)=[];     A(A>numel(i))=[];
    i(A)=1;
end

i=i((2+lr(1)):end-(1+lr(2)));



















