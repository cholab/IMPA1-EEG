function scalp(M,chanlocs_handle,drawmode,sc,dynamic,map,gridsc,LBL)
% scalp(M,chanlocs,drawmode,sc,dynamic,map,gridsc,LBL)
%To call direct plotting to specific handle:
% scalp(M,h,drawmode,sc,dynamic,map,gridsc)  [chanlocs needs to be global]


if ~nargin
    help scalp
else

if islogical(M)    
    M=double(M);
    sc=[.4 .6];
end
if size(M,1)==1
    M=M(:);
end

scmat=0;     
multiplot=0;
if isnumeric(chanlocs_handle) 
    global chanlocs
    tohandle=true;
    hh=chanlocs_handle;
else
    tohandle=false;    
    chanlocs=chanlocs_handle;
end
    
    
if ~exist('dynamic','var')||isempty(dynamic) 
    dynamic=0;
end
if numel(dynamic)>1
    ev=dynamic(2:end);    
    dynamic=dynamic(1);
else
    ev=[];
end
if dynamic<-1
    scmat=1;
end
if ~exist('gridsc','var') ||isempty(gridsc)
    gridsc=100;
end

if exist('drawmode','var') && isstruct(drawmode)
     AREA=drawmode;
     areamode=1;
     drawmode='areamode';
else
     areamode=0;
end


el=size(chanlocs,2); 
%interpret M
if isempty(M)
    showelectrodes=1;
    M=zeros(el,1);
elseif numel(M)<el && ~areamode
    showelectrodes=1;
    mt=M;
    M=zeros(el,1);
    M(mt(~isnan(mt)))=1;
    sc=[.4 .6];
else
    showelectrodes=0;
end
dim=ones(1,4); %up to 4d input
dim(1:numel(size(M)))=size(M);


eldim=find(el==dim);
if isempty(eldim)
    if dim(1)==dim(2)
       issurface=1; %it looks like a surface (defaults might be overriden never use 129x129 surfaces!)        
    else
        if ~areamode
        error('Input does not look like values from electrodes nor like a surface!')
        else
        issurface=0;  
        end
    end
else
issurface=0; 
if eldim(1)>1 && dim(4)==1 %electrodes are along d2-3 and it is not a 4d matrix
    scmat=1;         %this invokes scmatrix
end
end






if ~scmat
extrad=2+issurface;
if dim(extrad)==1 %no extra dimensions
multiplot=0;  dynamic=0; %n=t=1,single plot
else %One or Two extra dimensions...
if dim(extrad+1)==1 %One extra dimension
    if    dynamic; multiplot=0;   T=dim(extrad); %extra dim codes time
    else  multiplot=1;   nim=dim(extrad);  %extra dim codes multiple plots
    end
else %Two extra dimensions
if ~dynamic   
dynamic=.01; %overrides set value, 3D plotting is needed
end
multiplot=1;
T=dim(extrad); %first extra dim codes time
nim=dim(extrad+1);  %second extra dim codes multiple plots
end
end
end
if multiplot && nim>100 && nim>el && ~dynamic
%if extra dimension is very long and longer than the number of channels,
%multiplot is replaced by dynamic
        dynamic=.01; T=nim; multiplot=0;
end



% DRAWING MODES
mn=0; vn=0; projs=0; view=[];%matrix needed, vector needed
 if ~exist('drawmode','var')||isempty(drawmode) %setup defaults
     if issurface;      drawmode='contourf';   mn=1;
     else                  drawmode='scatter2d';  vn=1;    
     end
 else %intepret drawmodes...
    if ~areamode
     drawmode=lower(drawmode);
    if  ~isempty(strfind(drawmode,'3d'));  drawmode='scatter3d';   %there is ony a 3d plot so far  
    else
        if  ~isempty(strfind(drawmode,'pl'))  %combined plot is evoked
                  PLUS='PLUS';  vn=1; mn=1;
                  drawmode(strfind(drawmode,'pl'):end)=[];                  
        else   PLUS=[];
        end
        if  ~isempty(strfind(drawmode,'im'));      drawmode='imagesc'; mn=1;  
        elseif ~isempty(strfind(drawmode,'rf'));  drawmode='contourf';  mn=1;  
        elseif ~isempty(strfind(drawmode,'co')); drawmode='contour';  mn=1;  
        elseif ~isempty(strfind(drawmode,'pr')) 
            if ~isempty(strfind(drawmode,'all')) 
               projs=1; mn=1; vn=1;
           end
           view=drawmode(end);
           drawmode='scatterPROJ';                            
        else
                if       issurface; drawmode='contourf';  mn=1;
                else   drawmode='scatter2d'; vn=1;    
                end
        end
           drawmode=[drawmode PLUS];
    end
     end
 end
% ...and see if we have everything
if areamode %load specific locations for plotting(plocs) 
plocs=ch2plocs(chanlocs,AREA);     
scat_im=M; surf_im=[0 0 0 0];      
else
%load locations for plotting(plocs)    
plocs=ch2plocs(chanlocs,gridsc);    
plocs.view=view;
if projs; if ~dynamic; T=1; dynamic=.01; mn=1; end; end
if issurface
   surf_im=M; scat_im=[0 0 0];
   if vn;            scat_im=ch2surf(chanlocs,gridsc,surf_im);       
   end
else
   scat_im=M; surf_im=[0 0 0 0];      
   if mn;           surf_im=ch2surf(chanlocs,gridsc,scat_im);       
   end
end
end


% COLORMAPS 
 if ~exist('map','var')||isempty(map)
    map.im='jet';     map.scat='jet';    %defaults
 else
    tmap=map; clear map
    if iscell(tmap) && max(size(tmap)==2)  
        map.im=tmap{1}; map.scat=tmap{2};       
    else
        if iscell(tmap); tmap=tmap{1} ; end
        map.im=tmap;     map.scat=tmap;
    end
 end
 
%SCALES
if ~exist('sc','var')||isempty(sc)
    sc=quantile(M(~isnan(M)),[.05  .95]);
elseif numel(sc)==1 
    sc=[-abs(sc)  abs(sc)];
elseif numel(sc)==2
    sc=[sc(1) sc(2)]; 
end
v=linspace(sc(1),sc(2),100);
map.v=v([1 5 10 15   85 90 95 100]);
map.v=v([1:20:100 100]);
map.scale=sc;


%CALL PLOTTING FUNCTIONS
hold off; 
nm=1; vn=1; tm=1; tv=1;
if scmat
scmatrix([plocs.x plocs.y], M,sc)
elseif tohandle
nm=1; nv=1;    
for h=1:numel(hh)   
         if mn; nm=h; end;  if vn; nv=h; end   
         subplot(hh(h))        
         cla
         scalplot(plocs,scat_im(:,nv),surf_im(:,:,nm),drawmode,map,showelectrodes);    
end
else
if dynamic  
    if multiplot
    %MULTIPLE DYNAMIC PLOTS 
       switch drawmode 
         case 'scatter2d' %quick CData replacement
         POS=subplotter(nim);
         hh=subplotter(POS);
         for nn=1:nim
         POS.n=nn;
         hh(nn)=subplotter(POS);
         scalplot(plocs,scat_im(:,1,nn),surf_im,drawmode,map,showelectrodes);
         end
         lbwh=[min(POS.lm) min(POS.bm) max(POS.lm)+POS.w(find(max(POS.lm)==POS.lm,1,'first')) max(POS.bm)+POS.h(find(max(POS.bm)==POS.bm,1,'first'))];
         set(gcf,'color','k')
         if lbwh(2)==0; lbwh(2)=0.1; end
             
         annotation('Line', lbwh(1)+([0  lbwh(3)]), [lbwh(2)  lbwh(2) ] )
         vb=lbwh(4)/40; vb=lbwh(2)+[-vb vb];
         zb=lbwh(1); hb=lbwh(3);
         if ~isempty(ev)
           evt=zb+hb*ev;
           for nev=1:numel(ev)
           annotation('Line',[evt(nev)  evt(nev)], vb,'color','r','Linewidth',5)             
           end
         end
         cursor=annotation('Line',zb+[0 0], vb,'color','w','Linewidth',5);                     
         axis off
        
         for t=2:T;
         for nn=1:nim
         scalplot(hh(nn),scat_im(:,t,nn),map.scale,map.scat);
         set(cursor,'X',zb+hb*[t/T t/T])
         drawnow         
         end             

%          if dynamic>.05;  pause(dynamic);  end
         end
          
             otherwise %redraw (slow)
        POS=subplotter(nim);
        for t=1:T; if mn; tm=t; end ; if vn;  tv=t; end               
            for n=1:nim; if mn; nm=n; end;  if vn; nv=n; end                
             POS.n=n;   
    %          if t>1; cla(h(n)); end
             subplotter(POS);
             cla
             scalplot(plocs,scat_im(:,tv,nv),surf_im(:,tm,nm),drawmode,map,showelectrodes);
             drawnow         
             if dynamic>.05; pause(dynamic);  end
            end
        end        
        end
    else
    %SINGLE DYNAMIC PLOT        
     tA=0; tZ=1; tAZ=1;    
    if~isempty(ev)
        evt=repmat((tA+ev*tAZ)',[1 2]);        
    end
    switch drawmode 
     case 'scatter2d' %quick CData replacement
         scalplot(plocs,scat_im(:,1),surf_im,drawmode,map,showelectrodes);   
         hh=gca;   
         set(gcf,'color','k')
         lbwh=get(hh,'Position');         
         annotation('Line', lbwh(1)+([0  lbwh(3)]), [lbwh(2)  lbwh(2) ] )
         vb=lbwh(4)/40; vb=lbwh(2)+[-vb vb];
         zb=lbwh(1); hb=lbwh(3);
         if ~isempty(ev)
           evt=zb+hb*ev;
           for nev=1:numel(ev)
           annotation('Line',[evt(nev)  evt(nev)], vb,'color','r','Linewidth',5)             
           end
         end
         cursor=annotation('Line',zb+[0 0], vb,'color','w','Linewidth',5);                     
         axis off
        
         for t=2:T;              
         scalplot(hh,scat_im(:,t),map.scale,map.scat);
         set(cursor,'X',zb+hb*[t/T t/T])
         drawnow
         if dynamic>.05;  pause(dynamic);  end
         end
            
        otherwise %redrawing
            if projs
              POS=subplotter(9);
              POS=subplotter(POS);
              cla(POS(1)); cla(POS(3)); cla(POS(7)); cla(POS(9));
            end
    for t=1:T;   if mn; tm=t; end ; if vn;  tv=t; end      
%          if t>1; cla(gca); end
         if projs
             axes(POS(2)); plocs.view='a';
                 scalplot(plocs,scat_im(:,tv),surf_im(:,:,tm),drawmode,map,showelectrodes);
             axes(POS(8)); plocs.view='p';
                 scalplot(plocs,scat_im(:,tv),surf_im(:,:,tm),drawmode,map,showelectrodes);
             axes(POS(4)); plocs.view='l';
                 scalplot(plocs,scat_im(:,tv),surf_im(:,:,tm),drawmode,map,showelectrodes);
             axes(POS(6));  plocs.view='r';
                 scalplot(plocs,scat_im(:,tv),surf_im(:,:,tm),drawmode,map,showelectrodes);
             axes(POS(5));
                 scalplot(plocs,scat_im(:,tv),surf_im(:,:,tm),'contour',map,showelectrodes);             
         else
        
        if t==1
            h=gca;
        end
        subplot(h);
        scalplot(plocs,scat_im(:,tv),surf_im(:,:,tm),drawmode,map,showelectrodes);
        
        if t==1
         lbwh=get(gca,'Position');
         br=subplot('Position',[lbwh(1) lbwh(2)/2 lbwh(3) lbwh(2)/3]);
         axis off
        end
         subplot(br)
         line([0 1],[.5 .5],'LineStyle','-','color',[.5 .5 .5])
         line([t/T t/T],[0 1] ,'LineStyle','-','color','w','Linewidth',4)
         if ~isempty(ev)
            line(evt',[0 1] ,'LineStyle','-','color', 'r','Linewidth',4)          
         end
         end
         drawnow
         if dynamic>.05;  pause(dynamic);  end
    end
    end
    end
    
else %not dynamic
   if multiplot
    %MULTIPLE STATIC PLOTS       
    POS=subplotter(nim);
        for n=1:nim;  if mn; nm=n; end;  if vn; nv=n; end   
         POS.n=n;   
         subplotter(POS);   
         scalplot(plocs,scat_im(:,nv),surf_im(:,:,nm),drawmode,map,showelectrodes);
         if exist('LBL','var') && ~isempty(LBL)
             title(LBL{n})
         end
        end      
   else %n=t=1
         %SINGLE  STATIC PLOT
         scalplot(plocs,scat_im,surf_im,drawmode,map,showelectrodes);
   end
end
end
set(gcf,'color','k')
end


%    M(M<sc(2))=NaN;



            
function OUT=ch2surf(chanlocs,pl_gsc,in_im,surf_im)
%           surf_im=scalplot(chanlocs,gridsc,scat_im)
%           scat_im=scalplot(chanlocs,pl_gsc,surf_im)
if isstruct(pl_gsc) %scat_im -> surf_im
    d=mode(diff(pl_gsc.X));
    i=sub2ind(size(in_im),round(1+((pl_gsc.y-pl_gsc.Y(1))/d)) , round(1+((pl_gsc.x-pl_gsc.X(1))/d)));
    OUT=surf_im(i);
else %scat_im <- surf_im
    OUT=surfint(in_im,chanlocs,pl_gsc);   
end





