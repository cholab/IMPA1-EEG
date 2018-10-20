function [W,L,w,s,time,E]=decompose(X,dimensionality,method,binica_on,fastica_on,pca_on,cuda_on)
%     [W,L,w,s,time,E]=decompose(X,dimensionality,method,binica_on,fastica_on, pca_on,cuda_on)
%     [E,L]=decompose(X,'rank')
%Decompose X (channels,time,...) according to 'method'.
%Default is optimized call of ICA algorithms, [method=[]];
%       PCA(X) > L,PC
%       fastICA(PC) > w0 (initial guess)     [unless fastica_on==0, 1]
%       runICA*(PC) > w,s (final decomposition)  *[binICA if binica_on==1]
%Other accepted methods are pca rotated solutions ('varimax', promax'),
%anything else just calls pca.
%'dimensionality' reduce dimensionality to 'dimensionality%' variability or
%dimensionality components. If' 'rank' only explained variability and L are
%output, and the rank will be length(E) or size(L,2).
%* 'binica_on' calls binica, if it is a directory binica is called but files
%are stored in a temporary directory which is removed after computation.
%Recall that  activations will be:
%     ACT=W*X;             where W=w*s*pinv(L);
%And projections      
%     Pi=iW(:,i)*ACT(i,:)  where iW=pinv(W)


if ~nargin
    help decompose
    return
end


warning('off','MATLAB:RandStream:ActivatingLegacyGenerators')
warning('off','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState') 

if isa(X,'single')
    X=double(X);
end

justrank=0;
if exist('dimensionality','var') 
    if strcmpi(dimensionality,'rank')
    justrank=1;
    dimensionality=[];
    elseif dimensionality==0
        dimensionality=[];        
    end    
else
    dimensionality=[];
end

if ~isempty(dimensionality)
    pca_on=true;
end

binicadir=false;
if ~exist('binica_on','var') || isempty(binica_on)
    binica_on=false;
else
    if binica_on
    if ischar(binica_on)
    try [~,~]=mkdir(binica_on); end
    if isdir(binica_on)
        binicadir=binica_on;
        binica_on=true;
    else
       warning('Could not create a directory for binica') 
       binica_on=true;
    end
    end
    end
    binica_on=binica_on==1;
end
if binicadir
    cd(binicadir)
end

if ~exist('fastica_on','var') || isempty(fastica_on)
    fastica_on=true;
end
if ~exist('pca_on','var') || isempty(pca_on)
    pca_on=true;
end
if ~exist('cuda_on','var') || isempty(pca_on) %added by CPW 6-21-16
    cuda_on=false;
end

if ~exist('method','var') || isempty(method)
    ica_pca=1;
else
    ica_pca=2;
    if ~isempty(strfind(lower(method),'varimax'))
        rot='varimax';
    elseif  ~isempty(strfind(lower(method),'promax'))
        rot='promax';
    else
        rot=0;
    end
end



X=X(:,:);
nnans=~any(isnan(X));
display(['_' datestr(now) '_____COMPONENT ANALYSIS < decompose.m_____'])
tic; time=[];

if pca_on
%% PCA
fprintf('Going to principal components [PCA]... ') 
[L,PC,E]=princomp(X(:,nnans)'); time(1)=toc;
fprintf([' done.  [It took ' secs2str(time(1)) ']']) 
E=E/sum(E);

if exist('dimensionality','var') && ~isempty(dimensionality)
    fprintf('\nReducing dimensionality to')
    if ischar(dimensionality) %it is %of var explained
        dimensionality(dimensionality=='%')=[];
        dimensionality=str2double(dimensionality);
        ncomp=find(cumsum(E)>=(dimensionality/100),1,'first'); % **at least** this %var
        if isempty(ncomp); ncomp=0;end
         fprintf(' %s%% variability, i.e.',num2str(dimensionality))  
         if ncomp<2
          ncomp=2;   
          dimensionality=round(sum(E(1:2))*1000)/10;
         fprintf('...increasing  to  %s%% variability, i.e.',num2str(dimensionality))               
         end
    else
        ncomp=dimensionality;
    end
     fprintf(' the first %s components.', num2str(ncomp))           
else
    if fastica_on
    e=fastica(X(:,nnans),'verbose','off','only','pca');%see how many components fastica would retain    
    ncomp=size(e,2);
    fprintf('\nRetaining %s/%s components ',num2str(ncomp),num2str(size(L,2))) 
    else
      ncomp=size(PC,2);
    end
end
PC=PC(:,1:ncomp)';
L=L(:,1:ncomp); 
E=E(1:ncomp);    
end

if justrank
W=E; %first output is E--> length(E) is rank
fprintf(['\n']) 
else
if ica_pca==1
if fastica_on    
%% FASTICA
fprintf('\nLearning initial guess for component separation [FASTICA-tanh]... ')  
% [~,w0]=fastica(PC,'approach','symm','verbose', 'off','g','tanh','initGuess',wguess); 
if pca_on
[~,w0]=fastica(PC,'approach','symm','verbose', 'off','g','tanh'); %symm is faster, tanh handles better quasi gaussian components
else
[~,w0]=fastica(X(:,nnans),'approach','symm','verbose', 'off','g','tanh'); %symm is faster, tanh handles better quasi gaussian components    
end
time(2)=toc;
fprintf([' done.  [It took ' secs2str(time(2)-time(1)) ']']) 
else
    w0=false;
end

if ~binica_on
%% RUNICA
if pca_on
if w0
fprintf('\nRunning initialized extended ICA [RUNICA]... ')    
    if cuda_on  % added try-catch segment to attempt CUDAICA with data --CPW 6/21/16 
        try
            fprintf('\nAttempting CUDAICA first...')
            [w,s]=cudaica(PC,'weights',w0,'extended',1,'verbose','off');
        catch
            fprintf('\nCUDAICA failed, defaulting to RUNICA...')
            [w,s]=runica(PC,'weights',w0,'extended',1,'verbose','off');
        end
    else
        [w,s]=runica(PC,'weights',w0,'extended',1,'verbose','off');
    end
else % in case fastica did not converge
    if cuda_on  % added try-catch segment to attempt CUDAICA with data --CPW 6/21/16 
        try
            fprintf('\nAttempting CUDAICA first...')
            [w,s]=cudaica(PC,'extended',1,'verbose','off');
        catch
            fprintf('\nRunning NON initialized extended ICA [RUNICA]... ')        
            [w,s]=runica(PC,'extended',1,'verbose','off');   
        end
    else
        [w,s]=runica(PC,'extended',1,'verbose','off');   
    end
end
else
if w0
fprintf('\nRunning initialized extended ICA [RUNICA]... ')    
[w,s]=runica(X(:,nnans),'weights',w0,'extended',1,'verbose','off');
else % in case fastica did not converge
fprintf('\nRunning NON initialized extended ICA [RUNICA]... ')        
[w,s]=runica(X(:,nnans),'extended',1,'verbose','off');   
end    
end
% % % [~,s0]=fastica(PC,'verbose','off','only','white');
% % % [w2,s2]=runica(PC,'weights',s0\w0,'extended',1,'verbose','off');
else
if pca_on    
if w0
fprintf('\nRunning initialized extended ICA [BINICA]... ')    
[w,s]=binica(PC,'weights',w0,'extended',1,'verbose','off');
else % in case fastica did not converge
fprintf('\nRunning NON initialized extended ICA [BINICA]... ')        
[w,s]=binica(PC,'extended',1,'verbose','off');   
end    
else
if w0
fprintf('\nRunning initialized extended ICA [BINICA]... ')    
[w,s]=binica(X(:,nnans),'weights',w0,'extended',1,'verbose','off');
else % in case fastica did not converge
fprintf('\nRunning NON initialized extended ICA [BINICA]... ')        
[w,s]=binica(X(:,nnans),'extended',1,'verbose','off');   
end      
end
%cleanup intermediate files
% fprintf('\n deleting intermediate files... ')        
% delete binica*.sph    %put it into binica
% delete binica*.sc
% delete binica*.wts
% delete bias_after_adjust
if binicadir
fprintf('(deleting temporary directory)\n ')            
   binicadir=pwd;   
   cd(binicadir(1:find(binicadir=='/',1,'last')))
   [~,~]=rmdir(binicadir,'s');   
end
end
time(3)=toc;      
fprintf([' done.  [It took ' secs2str(time(3)-time(2)) ']\n']) 

%% BASIC ALGEBRA for data2components full transformation W:
if pca_on
W=w*s*pinv(L);
else
W=w*s;    
end
display(['_' datestr(now) '___________________________'])

else
   if rot
   fprintf('\nRotating...')    
   L=rotatefactors(L,'Method',rot,'normalize','off');  
   end
   W=pinv(L); s=nan(size(W)); w=s;
   fprintf([' done.\n']) 
end

if ica_pca==1 || ischar(rot)
try 
fprintf('Measuring explained variability and sorting components for oblique solutions... ') 
if pca_on
localE=1-sum(E);
else
localE=1;
end
A=data2act(X,W);
[~,E,i]=explainedvar(A,pinv(W));
E=E*localE;
E=E(i);
if pca_on
    L=L(:,i);
end
W=W(i,:);
s=s(:,i);
fprintf(' done.\n') 
catch
    save Wtemp_error W
end
end
end

  
 

 
% % CHECK DECOMPOSITION %this kind of works,  not in a sensible way though
% wlength=200;
% A=A(:,1:wlength*fix(size(A,2)/wlength));
% A=reshape(A,size(A,1),wlength,[]);
% bA=squeeze(binactivity(A,wlength));
% [~,i]=max(bA,[],2);
% i=unique(i);
% A=A(:,:,i);
% MI=hb_mi(A(:,:)');
% MI=transf(MI,'Z',2);
% [r,c]=find(MI>5);
% i=sort([r c],2);
% i=unique(i,'rows');
% ii=unique(i(:));
% for k=1:numel(ii)
%   seed=ii(k);    
%   [r c]=find(seed==i);
%   gr=i(r,:);
%   gr=unique(gr(:));  
%   i(r,:)=[];
%   GR(k,1:numel(gr))=gr;
% end
% GR(sum(GR,2)==0,:)=[];
% 
% W=zeros(size(w0));
% for n=1:size(GR,1)
%  gr=GR(n,:); gr(gr==0)=[];
%  W(n,:)=mean(w0(gr,:)); 
%  [~,im]=max(std(w0(gr,:),1,2));
%  W(n,:)=w0(gr(im),:);    
% end
% W2=w0(~ismember(1:size(A,1),ii),:);
% W(n+(1:size(W2,1)),:)=W2;




