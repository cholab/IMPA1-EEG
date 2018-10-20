function [S,Sz]=surfint(X,chanlocs,gridsc,progresslbl)
% S=surfint(X,chanlocs,gridsc)

if ~nargin 
    help surfint
else
if ~exist('gridsc','var')||isempty(gridsc)
    gridsc=67;
end
if isa(X,'single')
    X=double(X);
end



Th=deg2rad(cell2mat({chanlocs.theta})');
Rd=cell2mat({chanlocs.radius})';
[x,y] = pol2cart(Th,Rd);
% This to rotate nrad clockwise:
nrad=pi/2;
allcoords = (y + x*sqrt(-1))*exp(sqrt(-1)*nrad);
x = imag(allcoords);
y = real(allcoords);
el=numel(x);
xi = linspace(min(x),max(x),gridsc);   % x-axis description (row vector)
yi = linspace(min(y),max(y),gridsc);   % y-axis description (row vector)
[qx,qy] = meshgrid(xi,yi);

[n1,ns]=size(X);
if ns==el && n1~=el
    X=X';   [~,ns]=size(X);
end


if ns==1
     F = TriScatteredInterp(x,y,X);
     S = F(qx,qy); 
else
if ~exist('progresslbl','var')
    progresslbl=[];
end
S=zeros(gridsc,gridsc,ns);
% h = waitbar(0,[progresslbl 'surface rendering...']);
fprintf('%s',[progresslbl ' surface rendering '])
ns2=round(linspace(1,ns,10));
for n=ns2
    fprintf('%s', '.')   
end
for n=1:ns
if any(n==ns2)
    fprintf('\b\b%s', ' ')   
end
% [~,S(:,:,n)]=topoplot(X(:,n),chanlocs,'noplot','on');
% [Xi,Yi,S2(:,:,n)]= griddata(y,x,X(:,n),yi',xi,'invdist'); % interpolate data
% waitbar(n/ns,h)
F = TriScatteredInterp(x,y,X(:,n));
S(:,:,n) = F(qx,qy);
end
% fprintf('\b%c', 'done.')   
% display(' ')
% close(h)
fprintf('%s\n', ': done.')   
end
end

return

dim3=0;
x2=cell2mat({chanlocs.X})';
y2=cell2mat({chanlocs.Y})';
z2=cell2mat({chanlocs.Z})';

xi = linspace(min(x2)-1,max(x2)+1,gridsc);   % x-axis description (row vector)
yi = linspace(min(y2)-1,max(y2)+1,gridsc);   % y-axis description (row vector)
zi = linspace(min(y2),max(y2),gridsc);   % y-axis description (row vector)
[qx2,qy2] = meshgrid(xi,yi);
[~,qz2] = meshgrid(xi2,zi2);

F = TriScatteredInterp(x2,y2,z2,el);
qel2 = F(qx2,qy2,qz2); 


F = TriScatteredInterp(x2,y2,z2);
qz2 = F(qx2,qy2); 

[x,y,z] = sphere(50);
hgttilt = hgtransform;
hgrotate = hgtransform('parent',hgttilt);
% Set display properties
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'gouraud';
props.Cdata = topo;             % Use topo grid as a texturemap
props.Parent = hgrotate;              % Make hgtransform surface parent
hsurf = surface(x,y,z,props);   % Draw 3-D view
colormap(cmap)      % Use special terrain colormap defined above


[X,Y,Z] = meshgrid(xi2,yi2,zi2);
p = patch(isosurface(x,y,z,v,-3));


Th=deg2rad(cell2mat({chanlocs.theta})');
Rd=cell2mat({chanlocs.radius})';
[x,y] = pol2cart(Th,Rd);
xi = linspace(min(x),max(x),gridsc);   % x-axis description (row vector)

[faces,verts,colors] = isosurface(x,y,z,v,-3,x); 
patch('Vertices', verts, 'Faces', faces, ... 
    'FaceVertexCData', colors, ... 
    'FaceColor','interp', ... 
    'edgecolor', 'interp');
view(30,-15);
axis vis3d;
colormap copper[x,y,z,v] = flow; 
[faces,verts,colors] = isosurface(x,y,z,v,-3,x); 
patch('Vertices', verts, 'Faces', faces, ... 
    'FaceVertexCData', colors, ... 
    'FaceColor','interp', ... 
    'edgecolor', 'interp');
view(30,-15);
axis vis3d;
colormap copper[x,y,z,v] = flow; 
[faces,verts,colors] = isosurface(x,y,z,v,-3,x); 
patch('Vertices', verts, 'Faces', faces, ... 
    'FaceVertexCData', colors, ... 
    'FaceColor','interp', ... 
    'edgecolor', 'interp');
view(30,-15);
axis vis3d;
colormap copper










inspace(min(y),max(y),gridsc);   % y-axis description (row vector)
[qx,qy] = meshgrid(xi,yi);
F = TriScatteredInterp(x,y,X);
S = F(qx,qy); 

Fz = TriScatteredInterp(x2,y2,z2);
Sz = Fz(qx2,qy2); 
% surface(qx2,qy2,Sz,S,'FaceColor','texturemap','EdgeColor','none');
isosurface(qx2,qy2,Sz,S);
[x,y,z,v] = flow; 
[faces,verts,colors] = isosurface(X,Y,Z,v,-3,x); 
patch('Vertices', verts, 'Faces', faces, ... 
    'FaceVertexCData', colors, ... 
    'FaceColor','interp', ... 
    'edgecolor', 'interp');
view(30,-15);
axis vis3d;

% p=patch(surf2patch(qx2,qy2,Sz,S));
% shading faceted; view(3)
% set(p,'EdgeColor','none');









% Th=deg2rad(cell2mat({chanlocs.sph_theta})');
% Rd=cell2mat({chanlocs.radius})';
% Phi=deg2rad(cell2mat({chanlocs.sph_phi})');
% [x,y,z] = sph2cart(Th,Phi,Rd);
x=cell2mat({chanlocs.X})';
y=cell2mat({chanlocs.Y})';
z=cell2mat({chanlocs.Z})';

% This is to rotate nrad clockwise:
nrad=pi/2;
allcoords = (y + x*sqrt(-1))*exp(sqrt(-1)*nrad);
x = imag(allcoords);
y = real(allcoords);
el=numel(x);

xi = linspace(min(x),max(x),gridsc);   % x-axis description (row vector)
yi = linspace(min(y),max(y),gridsc);   % y-axis description (row vector)
zi = linspace(min(z),max(z),gridsc);   % z-axis description (row vector)

[qx,qy,qz] = meshgrid(xi,yi,zi);

[qx,qy] = meshgrid(xi,yi);
[qx,qz] = meshgrid(xi,zi);

[qx,qy,qz] = meshgrid(xi,yi,zi);
q = meshgrid(xi,yi,zi);


F = TriScatteredInterp(x,y,z,X);
S = F(qx,qy,qz); 
alpha('color')

load clown
surface(peaks,flipud(X),...
   'FaceColor','texturemap',...
   'EdgeColor','none',...
   'CDataMapping','direct')

colormap(map)
view(-35,45)












%   mask = (sqrt(Xi.^2 + Yi.^2)); % mask outside the plotting circle
%   mask=mask>=mask(1,round(gridsc/2));
%   mask=repmat(mask,[1 1 ns]);
%   S2(mask)=NaN;

  
  
%   
% tix = linspace(min(x),max(x),100);
% tiy = linspace(min(y),max(y),100);
% tiz = linspace(min(z),max(z),100);
% 
% [qx,qy] = meshgrid(xi,yi);
% qz = F(qx,qy);
% F = TriScatteredInterp([x y z],w);
% qw = F(qx,qy,qz);

  
  
  
  
  
  




% This is just to rotate:
%    allcoords = (y + x*sqrt(-1))*exp(sqrt(-1)*rotate);
%     x = imag(allcoords);
%     y = real(allcoords);

% gridsc = 67;
% xi = linspace(min(x),max(x),gridsc);   % x-axis description (row vector)
% yi = linspace(min(y),max(y),gridsc);   % y-axis description (row vector)
% 
%   [Xi,Yi,Zi] = griddata(y,x,intValues,yi',xi,'invdist'); % interpolate data
%   [Xi,Yi,Zi] = griddata(y,x,intContourVals,yi',xi,'invdist'); % interpolate data
%   
  
%   mask = (sqrt(Xi.^2 + Yi.^2) <= rmax); % mask outside the plotting circle
%   ii = find(mask == 0);
%   Zi(ii)  = NaN;                         % mask non-plotting voxels with NaNs  
%   ZiC(ii) = NaN; 
%   
%   
%   
% xyz=chanlocs;  
% xyz=[struct2mat(xyz,'X')  struct2mat(xyz,'Y') struct2mat(xyz,'Z')];
% xyz=round(xyz*100);
% x=xyz(:,1); y=xyz(:,2); z=xyz(:,3);
% 
% 
% 
% % Evaluate the interpolant at the locations (qx, qy). The corresponding value at these locations is qz .
% tix = linspace(min(x),max(x),100);
% tiy = linspace(min(y),max(y),100);
% tiz = linspace(min(z),max(z),100);
% 
% 
% 
% surface(qx,qy,qz,'FaceColor','texturemap', 'EdgeColor','none','FaceAlpha',.5)
% hold on;
% plot3(x,y,z,'o');
%  
% colormap jet
