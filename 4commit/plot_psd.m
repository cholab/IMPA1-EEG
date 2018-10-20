function [H]=plot_psd(f,PSD1,CI1,PSD2,CI2,AREA,log_on,savenames);

if ~exist('log_on','var')||isempty(log_on);
    log_on=false;
end
if ~exist('savenames','var')||isempty(savenames);
    saveflag=false;
else
    saveflag=true;
end

ticks=[find(round(f)==2,1,'first') find(round(f)==6,1,'first') find(round(f)==10,1,'first') find(round(f)==30,1,'first') find(round(f)==60,1,'first')];
tickLabel={'2'; '6'; '10'; '30'; '60'};

if log_on
   f=log10(f);
   PSD1=log10(PSD1);
   CI1=log10(CI1);
   PSD2=log10(PSD2);
   CI2=log10(CI2);
end


EB1=bsxfun(@minus,CI1,permute(PSD1,[1 2 4 3]));
EB1=EB1(:,:,[2 1],:); EB1(:,:,2,:)=abs(EB1(:,:,2,:)); 
EB2=bsxfun(@minus,CI2,permute(PSD2,[1 2 5 3 4]));
EB2=EB2(:,:,[2 1],:,:); EB2(:,:,2,:,:)=abs(EB2(:,:,2,:,:)); 

hi=max([CI1(:); CI2(:)]);
low=min([CI1(:); CI2(:)]);
if sign(low)==1 && sign(hi)==1;      % both positive
    ylims=[.9*low 1.1*hi];
elseif sign(low)==-1 && sign(hi)==1; % lower negative, upper positive
    
    ylims=[1.1*abs(low)*sign(low) 1.1*abs(hi)*sign(hi)];
else                                 % both negative
    ylims=[1.1*abs(low)*sign(low) .9*abs(hi)*sign(hi)];
end





dims=size(PSD2);
if numel(size(PSD2))>numel(size(PSD1))
    maxsubs=dims(end);
else maxsubs=1;
end
fn=1;
for cond=1:size(PSD1,3);
    for sub=1:maxsubs;
        figure('units','normalized','outerposition',[.1 .1 .8 .8]); 
        for j=1:8;

            if j<7; subplot(4,3,j)
            elseif j==7; subplot(4,3,8);
            else; subplot(4,3,11);
            end

            shadedErrorBar(f(3:end),PSD1(j,3:end,cond),squeeze(EB1(j,3:end,:,cond))',{'g','LineWidth',1.5,'markerfacecolor',[0 .8 0]});
            hold on;
            shadedErrorBar(f(3:end),PSD2(j,3:end,cond,sub),squeeze(EB2(j,3:end,:,cond,sub))',{'r','LineWidth',1.5,'markerfacecolor',[.8 0 0]});
            
            ylim(ylims);
            xlim([f(3) f(end)]);
            ax=gca;
            ax.XTick = f(ticks);
            ax.XTickLabel = tickLabel;
                       
            title(AREA(j).name);

        end
    
    if saveflag;
        export_fig(savenames{fn},'-r100')
        fn=fn+1;
    end

    end
end