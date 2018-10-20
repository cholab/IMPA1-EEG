function plotts(Y,yrange,i,focus,baseline,tl,mrks, dim,inlay)
%          plotts(Y,yrange,i,focus,baseline,tl,mrks, dim,inlay)
%   Y is 2D (if 3D, the third dimension is handled as trials)
%   yrange    sets vertical range, if 0 invokes zscore plotting [max(range(Y))]
%   i    are either channels (vector) or intervals ([onset offset]) to highlight
%   timeline   used for xtick  instead of indices
%   mrks       mark these points on timeline
if ~nargin
    help plotts 
    return
end
% it is a 2D plot, so 3D matrices are made 2D and third dimension is taken
% to encode trials
if ~exist('baseline','var') ||  isempty(baseline) 
    baseline=0;
end
if ~exist('focus','var') ||  isempty(focus)
    focus=0;
end
if baseline    
       Y=Y-repmat(mean(Y(:,baseline,:),2),[1 size(Y,2) 1]);
end
if ~exist('inlay','var') ||  isempty(inlay) 
    inlay=0;
end
ntr=size(Y,3); 
if ntr>1
    Y=Y(:,:);
    tr=linspace(0,size(Y,2),ntr+1);
end
% dim is time, usually the longest dimension
if ~exist('dim','var') || isempty(dim)
    dim=size(Y);
    [~,dim]=max(dim(1:2));
end

if dim==2 %time is expected to be the first dimension
    Y=Y';
end
[ntp,nch]=size(Y);

%if yrange is not provided everything should fit
if ~exist('yrange','var') || isempty(yrange)
   yrange=max(range(Y));
elseif yrange==0 %zero codes zscore plotting
    Y=zscore(Y);
    yrange=max(range(Y));
end
if yrange==0 
    yrange=1;
    display('flat lines?!')
end
if ~inlay
clf
colordef(gcf,'black')
set(gcf,'color',[.1 .1 .2])
colordef black
end

%this set channel offset based on vertical range for display
YY=Y+yrange*repmat(1:nch,ntp,1);

% h=axes;
axis([0 ntp 0  yrange*(nch+1)]);
line([0 ntp],repmat(yrange*(1:nch),[2 1] ),'LineStyle','-','color',.5*[.5 .5 1])
%line([0 ntp],repmat(yrange*(1:nch),[2 1] ),'LineStyle','-','color',[1 1 1])
hold on
%plot(YY,'w')
plot(YY,'defaultAxesColorOrder',colorcube)%use this for colors
xlim([0 ntp])
ylim([0 yrange*nch+yrange])
set(gca,'YTick',yrange*(1:nch),'Yticklabel',num2str((1:nch)'))

try
if focus==1 %default focus
        if ntr>5
            xlim([tr(1) tr(5+1)])
        elseif ntp>3000
            xlim([1 3000])    
        end      
elseif focus
      if ntr>1
          if focus<3
                 xlim([-tr(focus+2) tr(focus+1+3)])                            
          else
                 xlim([ tr(focus-3) tr(focus+1+3)])              
          end
      else
          if focus<3000
              xlim([-(focus+3000) focus+3000])
          else
              xlim([focus-3000  focus+3000])
          end
      end
end
end
pan xon

% 
% 
% uicontrol('Style', 'popup',...
%            'String', 'jet|hsv|hot|cool|gray',...
%            'Position', [],...
%            'Callback', @setmap);       % Popup function handle callback
%                                        % Implemented
%                                        
%                                        
%                                        
%                                        
%          rz=uicontrol('Style','pushbutton',...
%             'String','RAW/Z',...
%             'Callback','setrz');                                      
% %                                        


if exist('i','var') && ~isempty(i)
hold on
if isvector(i)    
    %vector of channels to higlight is provided
    if ~islogical(i) && numel(i)==nch && max(i)<2
        i=logical(i);
    end    
    plot(YY(:,i),'r')    
    
else  % interval to highlight is provided
%     if size(i,1)~=nch
%         i=i';
%     end
%     i(   i(:,1)<1  ,1)=1;
%     i(   i(:,2)>ntp  ,1)=ntp; 
%     for n=1:nch
%      plot(i(n,1):i(n,2),YY(i(n,1):i(n,2),n),'r')
%     end
    %here for 2 columns input=[onset offset]
    if size(i,1)<size(i,2)
        i=i';
    end
    for n=1:size(i,1)
     plot(i(n,1):i(n,2),YY(i(n,1):i(n,2),:),'r')
    end

end
end

if ntr>1
    line([tr(:) tr(:)], [0 yrange*(nch+1)], 'Color',.5*[0 0 1])
if exist('tl','var') && ~isempty(tl)
    it=round(linspace(1,ntp/ntr,5));
    it=repmat(it',[1 ntr])+repmat((0:ntr-1)*ntp/ntr,[5 1]);
    it=it(:); tl=repmat(tl,[1 ntr]);
    set(gca,'Xtick',it,'Xticklabel',tl(it))
end
else
if exist('tl','var') && ~isempty(tl)
if exist('mrks','var') && ~isempty(mrks)
%     [~,it]=ismember(mrks,tl);
    [~,it]=findind(mrks,tl);
    set(gca,'Xtick',it,'Xticklabel',tl(it))
    line([it; it], [0 yrange*(nch+1)], 'Color',.3*[0 1 0]) 
else    
    it=round(linspace(1,ntp,5));
    set(gca,'Xtick',it,'Xticklabel',tl(it))
end
end    
end


%display used range
title(num2str(yrange))
hold off










