function [H]=plotting_fx_MANUSCRIPT(prefix,frequencyband,grp)
%% This function quickly plots stuff for the IMPA1 manuscript
%   prefix is either 'raw' or 'whitened' to specify which analysis to load
%   frequencyband is the written-out band name with (e.g., 'theta')

colordef black

%% AMPA1 Plotting funtions for paper figures
% prefix=['raw']; %#ok<NASGU>
% prefix=['whitened']; %#ok<*NBRAK>
load([prefix '_qeeg.mat']);
[mMBP,~,~,seMBP]=groupav(AMeanPow,grp);
[mDFV,~,~,seDFV]=groupav(ADomFreqVar,grp);
xpos=[1:size(mMBP,1)]; %#ok<*NBRAK>
xshift=[-.2223 0 .2223];
regions={'Left Frontal'; 'Mid-Frontal'; 'Right Frontal'; ...
    'Left Temporal'; 'Central'; 'Right Temporal'; ...
    'Parietal'; 'Occipital'};
addpath('/home/chris/matlab/export_fig/')    
switch frequencyband
    case 'delta'
        B=1; name='Delta';
    case 'theta'
        B=2;   name='Theta';  
    case 'alpha'
        B=3;    name='Alpha'; 
    case 'beta'
        B=4;    name='Beta'; 
    case 'gamma1'
        B=5;     name='Gamma1';
    case 'gamma2'
        B=6;     name='Gamma2';    
end
%% Mean Band Power Results
% Plots show bar graphs of theta in each region
%
ylims=max(abs(stretch(mMBP(:,B,:,:),[2 1 3 4])))*1.3;

ffig; %figure; 
%eyes open
subplot(2,1,1); 
bar(squeeze(mMBP(:,B,1,:))); 
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mMBP(:,B,1,g)),squeeze(seMBP(:,B,1,g)),'Color','w','LineStyle','none','LineWidth',1.5);
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Power (uV^2)'); ylim([0 ylims]);  %Axis details
title([name ' Power -- Eyes Open']);
%eyes closed
subplot(2,1,2);
bar(squeeze(mMBP(:,B,2,:)));
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mMBP(:,B,2,g)),squeeze(seMBP(:,B,2,g)),'w','LineStyle','none','LineWidth',1.5);
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Power (uV^2)'); ylim([0 ylims]);  %Axis details
title([name ' Power -- Eyes Closed']);
export_fig(['/home/chris/Dropbox/IMPA1/fbat/' prefix '_' name '_Power_Fig_v1'],'-r200');

%% Dominant Frequency Variability 
ylims=max(abs(stretch(mDFV(:,B,:,:),[2 1 3 4])))*1.3;

ffig; %figure; 
%eyes open
subplot(2,1,1); 
bar(squeeze(mDFV(:,B,1,:))); 
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mDFV(:,B,1,g)),squeeze(seDFV(:,B,1,g)),'w','LineStyle','none','LineWidth',1.5);
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Mean Frequency Deviation (Hz)'); ylim([0 ylims]);  %Axis details
title([name ' Band Variability -- Eyes Open']);
%eyes closed
subplot(2,1,2);
bar(squeeze(mDFV(:,B,2,:)));
hold on; 
for g=1:size(unique(grp));
    errorbar(xpos+xshift(g),squeeze(mDFV(:,B,2,g)),squeeze(seDFV(:,B,2,g)),'w','LineStyle','none','LineWidth',1.5);
end
ax=gca; 
set(ax,'XTickLabel',regions)
ylabel('Mean Frequency Deviation (Hz)'); ylim([0 ylims]);  %Axis details
title([name ' Band Variability -- Eyes Closed']);
export_fig(['/home/chris/Dropbox/IMPA1/fbat/' prefix '_' name '_DFV_Fig_v1'],'-r200');
