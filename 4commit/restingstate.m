clear; clc; N=maxNumCompThreads(8);%max number of cores in grond
cd '/raid5/rcho/IMPA1/RAW/REST/';
datapath='/raid5/rcho/IMPA1/RAW/REST/';  
outpath='/raid5/rcho/IMPA1/MAT/REST/';  

addpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/'))
addpath(genpath('/raid5/rcho/TOOLS/EEGLAB_latest_vers/eeglab13_6_5b/'));
% addpath(genpath('/raid5/rcho/TOOLS/NIC/dpTOOLS/'));
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/fieldtrip-20150503/'));%removing overlapping functions
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/eeglab13_4_4b/functions/octavefunc/'));%removing overlapping functions
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/eeglab13_4_4b/plugins/Biosig3.0.7/'))
rmpath(genpath('/raid5/rcho/TOOLS/EEGLAB_latest_vers/'))
rmpath(genpath('/raid5/rcho/TOOLS/NIC/dpTOOLS/'))
rmpath(genpath('/raid5/rcho/TOOLS/NIC/MATLAB/eeglab12/'))
load chanlocsCOL
chanlocs=chanlocs([1:32 65]);
tmp=dir([datapath, '*.vhdr']);
fn={tmp(:).name}'; clear tmp

if length(fn)<30; 
    fn(26:30)=fn(25:29);
    fn{25}=[];
end

sr=5000; 
%t=-1:1/sr:5.3; %sampling rate and t is length of epoch in seconds 
for n=4:length(fn);
if ~isempty(fn{n})
    fn1=[outpath, 'impa', num2str(n,'%02i'), '_pre.mat'];
    
    if ~exist(fn1,'file');
        go=true;
    elseif ~ismember('W',who('-file',fn1)); 
        go=true; 
    else go=false;
    end

if go

    [X, ~, ~, Xstr]= loadraw2([datapath, fn{n}],sr);
    X=bsxfun(@minus, X, nanmean(X,2)); %detrending data, removal of overall mean from timeseries
    X=fqfiltering(X,sr,[],[59 61],150,0.2); %frequency filtering 

    [epochs1] =events2epoch('R100',Xstr); 
       if isempty(epochs1); [epochs1] =events2epoch('S 84',Xstr); end
       epochs1(:,1)=11;
    [epochs2] =events2epoch('R200',Xstr); 
       if isempty(epochs2); [epochs2] =events2epoch('S 64',Xstr); end
       epochs2(:,1)=12;
     
    template=zeros(48,2); template(:,2)=[0:12500:47*12500];
    epochs=[template; template] + [repmat(epochs1,48,1); repmat(epochs2,48,1)];
    epochs(:,3)=epochs(:,2)+12499;
    cond=epochs(:,1);%condition is first column 
    [v_raw2tr,v_tr2u,v_u2tr]=epochs2vect(epochs);

    ntr=size(epochs,1);ntp=1+unique(epochs(:,3)-epochs(:,2));

    X=reshape(X(:,v_raw2tr), [size(X,1), ntp, ntr]);
    [MASK,mask]=preicaPLUS(X, sr ,chanlocs,[],[],2);
    X=X(:,v_tr2u);
    
        save([fn1(1:end-4) '2.mat'], 'X','v_*','Xstr','epochs','cond','ntp','sr','MASK','mask');

    MASK=MASK(:,:,mask==4);
    btr=all(isnan(MASK));
    bch=all(isnan(MASK),2);
    X=ref(X);
    [~,~,~,u2TR,TR2u]=epochs2vect(epochs,[],~btr);
    X(find(bch),:,:)=[]; %#ok<FNDSB>

    X=X(:,u2TR);  
    %I need a u2D format    
    X=X(:,TR2u);  
    
    [W,L,w,s,time]=decompose(X,[],[],[],[],[],1);
    
        save(fn1,'-append','W','L','w','s','btr','bch');
        clear X epochs cond ntp MASK mask W L w s btr bch
        
end
end        
end 


return

tmp=dir([outpath, '*_pre.mat']);
fn={tmp(:).name}'; clear tmp
if length(fn)<30; 
    fn([26 27 29])=fn(25:27);
    fn{25}=[];
    fn{28}=[];
end



for n=1:length(fn);
if ~isempty(fn{n})
    fn1=[outpath, 'impa', num2str(n,'%02i'), '_pre'];
    fn2=[outpath, 'impa', num2str(n,'%02i'), '_post'];
    load(fn1,'X','epochs','MASK','sr','cond','mask','W')
    MASK=MASK(:,:,mask==4);
    bch=all(isnan(MASK),2);
    btr=all(isnan(MASK),1);
    [i,salvage]=posticaNEW(X,W,MASK,sr,epochs,chanlocs,[],[],1);
    find(i)'
    Xr=ref(X); Xr=reshape(Xr,[33 12500 96]); Xr=Xr(:,:,~btr); Xr=Xr(:,:);
    figure; plotts(Xr(~bch,:),150); % plot channel data
    figure; plotts(data2act(Xr(~bch,:),W'),10); % plot component activations
    figure; plotts(data2act(Xr(~bch,:),W(i,:)),10); % plot components identified by posticaNEW

    % Plot scalp topography with spectra
    figure; viscompFFT(W,cat(3,Xr(~bch,1:100000),Xr(~bch,end-99999:end)),chanlocs(~bch),[],[],'im',8,[],[],[],[],[1 sr]);
	% Plot scalp topography with time series
    figure; viscomp(W,cat(3,Xr(~bch,1:100000),Xr(~bch,end-99999:end)),chanlocs(~bch),[],[],'im',8);

    i2=find(i); i(:)=0; 
    i([i2' 4 26 18 8 29 28  ])=1;
    bad_ic=i;
        
    P=reproject(X,W,MASK,sr,epochs,chanlocs,bad_ic,3);
    figure; plotts(P(:,:),150);
    [X,btr,bch,MASK]=finaldatacleanup(P,sr,chanlocs,bch);
    TR=~btr;
    figure; plotts(X(:,:),150);
    save(fn2,'X','TR','cond','epochs','MASK','bch','bad_ic')
    restingPSDanalysis(X,sr,5,[5 6 9 13 19 20 24 25 29 ],cond(TR),1);
    
end
end