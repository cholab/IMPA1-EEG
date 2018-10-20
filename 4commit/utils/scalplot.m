function h=scalplot(plocs,in_im,surf_im,drawmode,map,showelectrodes)
%   h=scalplot(plocs,in_im,surf_im,drawmode,map,showelectrodes)
%   h=scalplot(h,scatvect,scale,map) [quick scatter2d CData replacement]

if isnumeric(plocs) %quick, plocs is an handle        
    if nargin<4              
        V=v2rgb(in_im,surf_im,'jet');         
    else
        V=v2rgb(in_im,surf_im,drawmode);                 
    end
    for n=1:numel(plocs)
    h=get(plocs(n),'Children');
    set(h(end), 'CData',V(:,:,n))  
    end
    
    
else
            if ~isempty(in_im) && size(in_im,2)~=3
                V=v2rgb(in_im,map.scale,map.scat);         
            end
            S=surf_im;          
            v=[];
            
         plus=strfind(drawmode,'PLUS') ;
         drawmode(plus:end)=[];
         d3=0;
        switch drawmode
            case 'scatter2d'
             if ~isempty(plus)
%                   colormap(loadcolormap(map.im))       
                  colormap(loadcolormap('bone'))                 
                  contour(plocs.X,plocs.Y,S,map.v)
                  hold on
             end   
            v=get(gca,'Position');v=v(3)*200;
            if v<40
            v=get(gca,'Position');v=v(3)*100;                
            end
            scatter(plocs.x,plocs.y,v,V,'filled')
            hold on
            
            if v>70 && ~showelectrodes
            scatter(plocs.x,plocs.y,v, 1*[1 1 1],'Linewidth',.01)
            end
            axis(gca,'off','square')     
            
            case 'scatter3d'                
            v=get(gca,'Position');v=v(3)*200;
            if v<40
            v=get(gca,'Position');v=v(3)*100;                
            end                
            v=v*5;
            scatter3(plocs.xyz(:,1),plocs.xyz(:,2),plocs.xyz(:,3),v,V,'filled')   
            hold on
            scatter3(plocs.xyz(:,1),plocs.xyz(:,2),plocs.xyz(:,3),v,.3*[1 1 1])                   
            axis(gca,'off','square')
            axis vis3d
            for dg=1:360
                view(dg,0);
                drawnow
            end
            
            d3=1;
            
            case 'scatterPROJ'
%          x post-ant, y is left-right
             pa=plocs.xyz(:,1);   
             lr=plocs.xyz(:,2);                
             y=plocs.xyz(:,3);     
              switch plocs.view
                  case 'a';    i=pa(:,1)>-.1;   x=lr; %frontal
                  case 'p' ;   i=pa(:,1)<.1;    x=lr; %posterior                       
                  case 'r' ;    i=lr(:,1)>-.1;    x=pa; %right side                                          
                  case 'l';     i=lr(:,1)<.1;     x=pa*-1; %left side                                            
              end              
              v=500;
            scatter(x(i),y(i),v, V(i,:),'filled')   
            
            hold on
            scatter(x(i),y(i),v, 0.3*[1 1 1])               
            axis(gca,'off','square')
            plocs.x=x(i); plocs.y=y(i);    
            case 'areamode'         
%             v=zscore(in_im); v=500*(.01+v-min(v));     
            v=1000;
            scatter(plocs.x,plocs.y,v,V,'filled')
            hold on
            scatter(plocs.x,plocs.y,v, .4*[1 1 1])
            axis(gca,'off','square')
                                        
            otherwise %matrix/surface images
                switch drawmode
            case 'imagesc'
            imagesc(plocs.X,plocs.Y,S,map.scale)
%             [x,y]=meshgrid(plocs.X,plocs.Y);
%             surfc(x,y,S)
            colormap(loadcolormap(map.im))           
            axis(gca,'xy')                   
            axis(gca,'off','square')
                 
           case 'contour'
           contour(plocs.X,plocs.Y,S,map.v)
           colormap(loadcolormap(map.im))
           axis(gca,'off','square')

           case 'contourf'
           contourf(plocs.X,plocs.Y,S,map.v)
           colormap(loadcolormap(map.im))
           axis(gca,'off','square')
           
           
                end
           if ~isempty(plus)
               hold on
               scatter(plocs.x,plocs.y,[],V,'filled')
               scatter(plocs.x,plocs.y,[],repmat(.4,[1 3]))
           end           
        end


if showelectrodes
if d3
    try
text(plocs.xyz(:,1),plocs.xyz(:,2),plocs.xyz(:,3),plocs.lbl,'HorizontalAlignment','center')        
    catch err
        display('no readable electrodes labels')
        display(err.message)        
text(plocs.xyz(:,1),plocs.xyz(:,2),plocs.xyz(:,3),num2str((1:numel(plocs.x))'),'HorizontalAlignment','center')
    end
else
    try
text(plocs.x,plocs.y,plocs.lbl,'HorizontalAlignment','center')                
    catch err
        display('no readable electrodes labels')        
        display(err.message)                
text(plocs.x,plocs.y,num2str((1:numel(plocs.x))'),'HorizontalAlignment','center')        
    end
end
end
end
