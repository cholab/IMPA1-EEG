%% IMPA1 Plotting funtions for paper figures
prefix=['raw']; %#ok<NASGU>
%prefix=['whitened']; %#ok<*NBRAK>
load([prefix '_qeeg.mat']);
[mMBP,~,~,seMBP]=groupav(AMeanPow,grp);
[mDFV,~,~,seDFV]=groupav(ADomFreqVar,grp);
xpos=[1:size(mMBP,1)];
xshift=[-.2223 0 .2223];
regions={'Left Frontal'; 'Mid-Frontal'; 'Right Frontal'; ...
    'Left Temporal'; 'Central'; 'Right Temporal'; ...
    'Parietal'; 'Occipital'};
addpath('/home/chris/matlab/export_fig/')    

addpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/'));
addpath(genpath('/raid5/rcho/TOOLS/CHRIS/'));
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/fieldtrip-20150503/'));%removing overlapping functions
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/eeglab13_4_4b/functions/octavefunc/'));%removing overlapping functions
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/eeglab13_4_4b/plugins/Biosig3.0.7/'));
rmpath(genpath('/raid5/rcho/TOOLS/EEGLAB_latest_vers/'));
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/eeglab12/'));

%% THETA POWER - MR edit

% c=1 eyes open
% c=2 eyes closed

for c=1:2
    ffig;
    subplot(20,3,[7 10 13 16]);
    MBPtest=squeeze(mMBP(1,2,c,:));
    hb=bar(MBPtest);
    hold on
    bar(1,MBPtest(1),'facecolor',[.08 .17 .55]);
    bar(2,MBPtest(2),'facecolor',[0 .75 .75]);
    bar(3,MBPtest(3),'facecolor',[1 1 0]);
    ylim([0 0.035]); ylabel('Power (uV^2)');
    xlabel('Left Frontal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mMBP(1,2,c,g)),squeeze(seMBP(1,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    subplot(20,3,[9 12 15 18]);
    MBPtest=squeeze(mMBP(3,2,c,:));
    hb=bar(MBPtest);
    hold on
    bar(1,MBPtest(1),'facecolor',[.08 .17 .55]);
    bar(2,MBPtest(2),'facecolor',[0 .75 .75]);
    bar(3,MBPtest(3),'facecolor',[1 1 0]);
    ylim([0 0.035]);
    xlabel('Right Frontal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mMBP(3,2,c,g)),squeeze(seMBP(3,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    for p=[2,5,8,11]
        subplot(4,3,p);
        if p==2
            sp=2;
        end
        if p==5
            sp=5;
        end
        if p==8
            sp=7;
        end
        if p==11
            sp=8;
        end
        MBPtest=squeeze(mMBP(sp,2,c,:));
        hb=bar(MBPtest);
        hold on
        bar(1,MBPtest(1),'facecolor',[.08 .17 .55]);
        bar(2,MBPtest(2),'facecolor',[0 .75 .75]);
        bar(3,MBPtest(3),'facecolor',[1 1 0]);
        ylim([0 0.035]);
        set(gca,'xticklabel',{[]})
        if p==2
            xlabel('Mid-Frontal');
        end
        if p==5
            xlabel('Central');
        end
        if p==8
            xlabel('Parietal');
        end
        if p==11
            xlabel('Occipital');
            ylabel('Power (uV^2)');
        end
        hold on;
        for g=1:size(unique(grp))
            errorbar(g,squeeze(mMBP(sp,2,c,g)),squeeze(seMBP(sp,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
        end
    end
    subplot(20,3,[25 28 31 34]);
    MBPtest=squeeze(mMBP(4,2,c,:));
    hb=bar(MBPtest);
    hold on
    bar(1,MBPtest(1),'facecolor',[.08 .17 .55]);
    bar(2,MBPtest(2),'facecolor',[0 .75 .75]);
    bar(3,MBPtest(3),'facecolor',[1 1 0]);
    ylim([0 0.035]); ylabel('Power (uV^2)');
    xlabel('Left Temporal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mMBP(4,2,c,g)),squeeze(seMBP(4,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    subplot(20,3,[27 30 33 36]);
    MBPtest=squeeze(mMBP(6,2,c,:));
    hb=bar(MBPtest);
    hold on
    bar(1,MBPtest(1),'facecolor',[.08 .17 .55]);
    bar(2,MBPtest(2),'facecolor',[0 .75 .75]);
    bar(3,MBPtest(3),'facecolor',[1 1 0]);
    ylim([0 0.035]);
    xlabel('Right Temporal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mMBP(6,2,c,g)),squeeze(seMBP(6,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    %fig=gcf;
    %set(fig,'name','Theta Power -- Eyes Open');
    %set(gcf,'color',[.9 .9 .9]);
end


%export_fig(['/home/chris/Dropbox/IMPA1/fbat/' prefix 'Theta_Power_Fig_v1'],'-r200');


%% DOMINANT FREQUENCY VARIABILITY - THETA MR EDIT

for c=1:2
    ffig;
    subplot(20,3,[7 10 13 16]);
    DFVtest=squeeze(mDFV(1,2,c,:));
    hb=bar(DFVtest);
    hold on
    bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
    bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
    bar(3,DFVtest(3),'facecolor',[1 1 0]);
    ylim([0 1.5]); %ylabel('Mean Frequency Deviation (Hz)');
    xlabel('Left Frontal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mDFV(1,2,c,g)),squeeze(seDFV(1,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    subplot(20,3,[9 12 15 18]);
    DFVtest=squeeze(mDFV(3,2,c,:));
    hb=bar(DFVtest);
    hold on
    bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
    bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
    bar(3,DFVtest(3),'facecolor',[1 1 0]);
    ylim([0 1.5]);
    xlabel('Right Frontal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mDFV(3,2,c,g)),squeeze(seDFV(3,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    for p=[2,5,8,11]
        subplot(4,3,p);
        if p==2
            sp=2;
        end
        if p==5
            sp=5;
        end
        if p==8
            sp=7;
        end
        if p==11
            sp=8;
        end
        DFVtest=squeeze(mDFV(sp,2,c,:));
        hb=bar(DFVtest);
        hold on
        bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
        bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
        bar(3,DFVtest(3),'facecolor',[1 1 0]);
        ylim([0 1.5]);
        set(gca,'xticklabel',{[]})
        if p==2
            xlabel('Mid-Frontal');
        end
        if p==5
            xlabel('Central');
        end
        if p==8
            xlabel('Parietal');
        end
        if p==11
            xlabel('Occipital');
            %ylabel('Mean Frequency Deviation (Hz)');
        end
        hold on;
        for g=1:size(unique(grp))
            errorbar(g,squeeze(mDFV(sp,2,c,g)),squeeze(seDFV(sp,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
        end
    end
    subplot(20,3,[25 28 31 34]);
    DFVtest=squeeze(mDFV(4,2,c,:));
    hb=bar(DFVtest);
    hold on
    bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
    bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
    bar(3,DFVtest(3),'facecolor',[1 1 0]);
    ylim([0 1.5]); ylabel('Mean Frequency Deviation (Hz)');
    xlabel('Left Temporal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mDFV(4,2,c,g)),squeeze(seDFV(4,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    subplot(20,3,[27 30 33 36]);
    DFVtest=squeeze(mDFV(6,2,c,:));
    hb=bar(DFVtest);
    hold on
    bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
    bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
    bar(3,DFVtest(3),'facecolor',[1 1 0]);
    ylim([0 1.5]);
    xlabel('Right Temporal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mDFV(6,2,c,g)),squeeze(seDFV(6,2,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    %fig=gcf;
    %set(fig,'name','Theta Variability -- Eyes Open');
    %set(gcf,'color',[1 1 1]);
end



%% DOMINANT FREQUENCY VARIABILITY - ALPHA MR EDIT

for c=1:2
    ffig;
    subplot(20,3,[7 10 13 16]);
    DFVtest=squeeze(mDFV(1,3,c,:));
    hb=bar(DFVtest);
    hold on
    bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
    bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
    bar(3,DFVtest(3),'facecolor',[1 1 0]);
    ylim([0 1.5]); %ylabel('Mean Frequency Deviation (Hz)');
    xlabel('Left Frontal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mDFV(1,3,c,g)),squeeze(seDFV(1,3,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    subplot(20,3,[9 12 15 18]);
    DFVtest=squeeze(mDFV(3,3,c,:));
    hb=bar(DFVtest);
    hold on
    bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
    bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
    bar(3,DFVtest(3),'facecolor',[1 1 0]);
    ylim([0 1.5]);
    xlabel('Right Frontal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mDFV(3,3,c,g)),squeeze(seDFV(3,3,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    for p=[2,5,8,11]
        subplot(4,3,p);
        if p==2
            sp=2;
        end
        if p==5
            sp=5;
        end
        if p==8
            sp=7;
        end
        if p==11
            sp=8;
        end
        DFVtest=squeeze(mDFV(sp,3,c,:));
        hb=bar(DFVtest);
        hold on
        bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
        bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
        bar(3,DFVtest(3),'facecolor',[1 1 0]);
        ylim([0 1.5]);
        set(gca,'xticklabel',{[]})
        if p==2
            xlabel('Mid-Frontal');
        end
        if p==5
            xlabel('Central');
        end
        if p==8
            xlabel('Parietal');
        end
        if p==11
            xlabel('Occipital');
            %ylabel('Mean Frequency Deviation (Hz)');
        end
        hold on;
        for g=1:size(unique(grp))
            errorbar(g,squeeze(mDFV(sp,3,c,g)),squeeze(seDFV(sp,3,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
        end
    end
    subplot(20,3,[25 28 31 34]);
    DFVtest=squeeze(mDFV(4,3,c,:));
    hb=bar(DFVtest);
    hold on
    bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
    bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
    bar(3,DFVtest(3),'facecolor',[1 1 0]);
    ylim([0 1.5]); ylabel('Mean Frequency Deviation (Hz)');
    xlabel('Left Temporal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mDFV(4,3,c,g)),squeeze(seDFV(4,3,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    subplot(20,3,[27 30 33 36]);
    DFVtest=squeeze(mDFV(6,3,c,:));
    hb=bar(DFVtest);
    hold on
    bar(1,DFVtest(1),'facecolor',[.08 .17 .55]);
    bar(2,DFVtest(2),'facecolor',[0 .75 .75]);
    bar(3,DFVtest(3),'facecolor',[1 1 0]);
    ylim([0 1.5]);
    xlabel('Right Temporal');
    set(gca,'xticklabel',{[]})
    hold on;
    for g=1:size(unique(grp))
        errorbar(g,squeeze(mDFV(6,3,c,g)),squeeze(seDFV(6,3,c,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
    end
    fig=gcf;
    set(fig,'name','Theta Variability -- Eyes Open');
    set(gcf,'color',[1 1 1]);
end




%% Mean Band Power Results -- Theta 
% Plots show bar graphs of theta in each region

ffig; %figure; 
%eyes open
subplot(2,1,1); 
bar(squeeze(mMBP(:,2,1,:))); 
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mMBP(:,2,1,g)),squeeze(seMBP(:,2,1,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Power (uV^2)'); ylim([0 0.04]);  %Axis details
title('Theta Power -- Eyes Open');
%eyes closed
subplot(2,1,2);
bar(squeeze(mMBP(:,2,2,:)));
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mMBP(:,2,2,g)),squeeze(seMBP(:,2,2,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Power (uV^2)'); ylim([0 0.04]);  %Axis details
title('Theta Power -- Eyes Closed');
%export_fig(['/home/chris/Dropbox/IMPA1/fbat/' prefix 'Theta_Power_Fig_v1'],'-r200');


%% Dominant Frequency Variability -- Theta

ffig; %figure; 
%eyes open
subplot(2,1,1); 
bar(squeeze(mDFV(:,2,1,:))); 
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mDFV(:,2,1,g)),squeeze(seDFV(:,2,1,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Mean Frequency Deviation (Hz)'); ylim([0 1.5]);  %Axis details
title('Theta Band Variability -- Eyes Open');
%eyes closed
ffig;
subplot(2,1,2);
bar(squeeze(mDFV(:,2,2,:)));
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mDFV(:,2,2,g)),squeeze(seDFV(:,2,2,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Mean Frequency Deviation (Hz)'); ylim([0 1.5]);  %Axis details
title('Theta Band Variability -- Eyes Closed');
%export_fig(['/home/chris/Dropbox/IMPA1/fbat/' prefix 'Theta_DFV_Fig_v1'],'-r200');


%% Dominant Frequency Variability -- Alpha

ffig; %figure; 
%eyes open
subplot(2,1,1); 
bar(squeeze(mDFV(:,3,1,:))); 
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mDFV(:,3,1,g)),squeeze(seDFV(:,3,1,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Mean Frequency Deviation (Hz)'); ylim([0 1.5]);  %Axis details
title('Alpha Band Variability -- Eyes Open');
%eyes closed
subplot(2,1,2);
bar(squeeze(mDFV(:,3,2,:)));
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mDFV(:,3,2,g)),squeeze(seDFV(:,3,2,g)),'LineStyle','none','LineWidth',1.5,'Color','k');
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Mean Frequency Deviation (Hz)'); ylim([0 1.5]);  %Axis details
title('Alpha Band Variability -- Eyes Closed');
export_fig(['/home/chris/Dropbox/IMPA1/fbat/'  prefix 'Alpha_DFV_Fig_v1'],'-r200');

