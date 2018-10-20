function [OUT,LBL]=el2a(IN,AREA,dim, modeswitch)
%From electrodes to areas (or from areas to electrodes), i.e. extract
%areawise averages (or map areas in electrode space):
%       OUT=el2a(IN,AREA,dim,modeswitch)
%modeswitch can be [0]:
%       0 el2a, averages
%       1 el2a, modeswitch 
%      -1 a2el, 129 electrodes
%      -n a2el, n electrodes

if ~nargin
    help  el2a
    return
end
   
    
if ~exist('modeswitch','var')||isempty(modeswitch)
modeswitch=0;
end

if ~exist('AREA','var')||isempty(AREA)
help el2a
error('Missing AREA information')
end

if isstruct(AREA)
    if nargout==2
    LBL={AREA.name}';
    end
    AREA={AREA.electrodes};    
end
AREA=AREA(:);
na=size(AREA,1);  

if isempty(IN)    
    modeswitch=-1;
    IN=(1:na)';
end

dims=size(IN);

if modeswitch>-1 %from electrodes 2 areas
   
if ~exist('dim','var')||isempty(dim)
  dim=find(dims==129);
if isempty(dim)
    dim=find(dims~=1);
    if numel(dim)>1    
     [~,dim]=min(dims);
     end
end
end

 

if dim==3
if numel(dims)==3
    dims=[dims 1];
end
OUT=zeros([dims(1:2) na dims(4)]);  
for d4=1:dims(4)
for a=1:na    
   OUT(:,:,a,d4)=averager(IN(:,:,AREA{a},d4),3,modeswitch);
end
end



elseif dim == 4
OUT=zeros([dims(1:3) na]);  
for d3=1:dims(3)
for a=1:na    
   OUT(:,:,d3,a)=averager(IN(:,:,d3,AREA{a}),4,modeswitch);
end
end
    

else
if numel(dims)<3
    dims(3)=1;
end
if numel(dims)<4
    dims(4)=1;
end

if dim==1
OUT=zeros([na dims(2:4)]);  
for nd4=1:dims(4)
for nd3=1:dims(3)
for a=1:na    
   OUT(a,:,nd3,nd4)=averager(IN(AREA{a},:,nd3,nd4),1,modeswitch);
end     
end
end

elseif dim==2
if numel(dims)<3
    dims(3)=1;
end
if numel(dims)<4
    dims(4)=1;
end
    
OUT=zeros([dims(1) na dims(3:4)]);  
for nd4=1:dims(4)
for nd3=1:dims(3)
for a=1:na    
   OUT(:,a,nd3)=averager(IN(:,AREA{a},nd3,nd4),2,modeswitch);
end     
end
end

end 

end


else %from areas 2 electrodes

    
if ~exist('dim','var')||isempty(dim)
  dim=find(dims==na);
  dim=dim(1);
end
if modeswitch==-1
  el=129;  
else
  el=abs(modeswitch);  
end
        
switch dim
 case 1
    OUT=nan([el dims(2:end)]);
    for n=1:na
        OUT(AREA{n},:,:,:)=repmat(IN(n,:,:,:),[numel(AREA{n}) 1 1 1]);        
    end
 case 2
    OUT=nan([dims(1) el dims(3:end)]);
    for n=1:na
        OUT(:,AREA{n},:,:)=repmat(IN(:,n,:,:),[1 numel(AREA{n})  1 1]);        
    end  
end

end
