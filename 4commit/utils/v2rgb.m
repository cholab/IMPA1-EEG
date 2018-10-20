function rgb=v2rgb(V,sc,map)
% rgb=v2rgb(V,sc,map)

if ~nargin
    help v2rgb
else
[nv np]=size(V);
if ~exist('sc','var')||isempty(sc)
i=V;    
mn=min(V); mx=max(V);    
for n=1:np
r=linspace(mn(n),mx(n),64);
[~,i(:,n)]=histc(V(:,n), r);
end
else
    if numel(sc)==1 
        if sc==1
             v=V(:);
             sc=quantile(v(~isnan(v)),[.1 .9]);
        else
            sc=[-abs(sc) abs(sc)];
        end
    end
    r=linspace(sc(1),sc(2),64);
    V(V<sc(1))=sc(1);    V(V>sc(2))=sc(2);    
     [~,i]=histc(V, r);
end
if ~exist('map','var') || isempty(map)
    map='Jet';
end
if ischar(map);
    map=loadcolormap(map);
end
nn=(i==0);
i(nn)=1;
rgb=zeros(nv,3,np);
rgb(:,1,:)=reshape(map(i(:),1),[nv,1,np]);
rgb(:,2,:)=reshape(map(i(:),2),[nv,1,np]);
rgb(:,3,:)=reshape(map(i(:),3),[nv,1,np]);
nn=repmat(reshape(nn,[nv,1,np]),[1 3 1]);
rgb(nn)=NaN;
end


