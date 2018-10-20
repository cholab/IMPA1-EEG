function OUT=subplotter(IN,space_ax,covered_width,covered_height,makeas1)
%To create the structure POS to place an array of npf subplots:
%   POS=subplotter(npf,subspace,covered_width,covered_height,makeas1)
%To creates the structure POS to center subplots at (X,Y):
%   POS=subplotter([X Y],subspace,[w h])
%To plot at POS and get the related handles:
%   h=subplotter(POS,axisoff)

if ~nargin
    help subplotter
    return
end
    
  %intepret IN to call function in correct sbmode    
  if isstruct(IN)%sublot from an existent POS structure 
     sbmode=2;
  else %create  a POS structure that positions subplots in a figure
       switch numel(IN)
           case  1 %nfg
               sbmode=1;
           case 2 
               sbmode=1;
               r=IN(1); c=IN(2);
               otherwise
               sbmode=3; 
               X=IN(:,1); Y=IN(:,2);
       end
  end
       
if sbmode==2 % sublot from an existent POS structure
 POS=IN;
 figure(gcf)
 if isfield(POS,'n') && ~isempty(POS.n)  %just one
 nl=POS.n-(fix(POS.n/POS.maximsperfig)*POS.maximsperfig);
 nl(nl==0)=POS.maximsperfig;
%  if nl==1 && POS.n>1
%  figure
%  end
%  subplot(POS.r,POS.c,nl)   
 OUT=subplot('Position',[POS.lm(nl) POS.bm(nl) POS.w(nl) POS.h(nl)]);
 box off
 else
 OUT=zeros(1,POS.npf); 
 for n=1:POS.npf     
 nl=n-(fix(n/POS.maximsperfig)*POS.maximsperfig);
 nl(nl==0)=POS.maximsperfig;
 if n>1 && nl==1 && POS.npf>1
 figure
 end
%  subplot(POS.r,POS.c,nl)   
 OUT(n)=subplot('Position',[POS.lm(n) POS.bm(n) POS.w(n) POS.h(n)]);
 box off
 if ~exist('space_ax','var')  || isempty(space_ax) || space_ax
     axis([0 1  0 1])
     text(.48,.48,num2str(n),'HorizontalAlignment','center','color','k')
     text(.52,.52,num2str(n),'HorizontalAlignment','center','color','w')     
     axis off
 end
 end      
end  
    
else    
  
    
    
    
    
%create  a POS structure that positions subplots in a figure
if sbmode==1 %just npf is provided
npf=IN;
POS.npf=npf;

maximsperfig=1000; %??
if npf>maximsperfig
    display(['Displaying only ' num2str(maximsperfig) '/' num2str(npf) ' axes!'])
end
npf(npf>maximsperfig)=maximsperfig;   
if ~exist('r','var')
r=floor(sqrt(npf));
c=ceil(npf/r);
end
maximsperfig=r*c;
POS.r=r;
POS.c=c;

POS.maximsperfig=maximsperfig;
%set axis dimensions/positions (they will cover a part (covered_width/height) of [0 1])
if ~exist('covered_width','var')||isempty(covered_width) 
     covered_width=1;
end
if ~exist('covered_height','var')||isempty(covered_height) 
     covered_height=covered_width;
end
w=covered_width/c; %width
h=covered_height/r; %height
% left margin and bottom margins
minlm=(1-covered_width)/(c+1);
lms=minlm+(0:c-1)*(minlm+w);
minbm=(1-covered_height)/(r+1);
bms=minbm+((r-1):-1:0)*(minbm+h);
lm=repmat(lms,[1 r])';
bm=repmat(bms',[1 c])';
bm=bm(:);
elseif sbmode==3  %subplots centered at (X,Y)
 npf=numel(X) ; 
 bestmosaic=0;
 if exist('covered_width','var') && ~isempty(covered_width) 
 wh=covered_width;     %here covered_width is actually [w h] 
     if numel(wh)==1 %it is w/h, but since it is redundant... 
         if wh==1 %...wh=1 invokes bestmosaic mode...
             bestmosaic=1; %...which forces alignment to reduce blank space
         else %otherwise it is a ratio of horizontal to vertical dimension
             rt=covered_width;
         end
     else
    w=wh(1);  h=wh(2); 
     end
 else %try to guess appropriate w and h
     w=min(pdist([X Y]))/1.4143; %for square inlays
     h=w;
 end
if bestmosaic
%  overlay a grid and see when I have no overlapping plots
r=numel(unique(round(10*(X/range(X))))); %initial guess
c=numel(unique(round(10*(Y/range(Y)))));    
nu=0;
while nu<npf; %while different nodes are mapped onto the same next coordinate
r=r+1; c=c+1;
x=linspace(min(X),max(X),c); 
y=linspace(min(Y),max(Y),r);
xy=[repmat(x,[1 r])'  reshape(repmat(y,[c 1]),r*c,1) ]; 
g=size(xy,1); 
dist=squareform(pdist([[X Y]; xy]));
[~,i]=min(dist(1:npf,end-g+1:end),[],2);
nu=numel(unique(i));
end
%can it be squeezed? (at the cost of altering relative positions)
x=1:c;  y=1:r;
xy=[repmat(x,[1 r])'  reshape(repmat(y,[c 1]),r*c,1) ]; 
f=zeros(r,c);  f(sub2ind(size(f), xy(i,2),xy(i,1)))=1;
f(sub2ind(size(f), xy(i,2),xy(i,1)))=1:npf; %locate original nodes in the new grid!
for nr=1:r-1
    i0=find(f(nr,:));     i1=find(f(nr+1,:));
if isempty(intersect(i0,i1))
    f(nr+1,i0)= f(nr,i0);
    f(nr,:)=0;
end
end
f=f(any(f,2),:);
for nc=1:c-1
    i0=find(f(:,nc));     i1=find(f(:,nc+1));
if isempty(intersect(i0,i1))
    f(i0,nc+1)= f(i0,nc);
    f(:,nc)=0;
end
end
f=f(:,any(f));
[~,c]=size(f);
c2=max(sum(f>0,1)); %this the best minimum n of columns to fit all the plots as adjacent plots
if c2<c-1
r2c=find(sum(f>0,2)==c2 | f(:,1)>0);
for nr=1:numel(r2c)
   j=r2c(nr);
   lf=f(j,:);
   i=find(f(j,:)==0);
   ni=fix(numel(i)/2);
   lf(i([1:ni end-ni+1:end]))=[];
   f(j,:)=0;
   f(j, ((c-numel(lf))/2 +(1:numel(lf))))=lf;
end
end
f=f(:,any(f));
[r,c]=size(f);

% if it is still suboptimal (because of symmetry requirement), let's make  positions independent (i.e. break the matrix)
c2=max(sum(f>0,1));
r2=max(sum(f>0,2));
if r2<r
if r-r2>c-c2 %if we gain more compressing rows
xb=zeros(sum(f(:)>0),1);
yb=zeros(sum(f(:)>0),1);
iif=zeros(sum(f(:)>0),1);
k1=0;
for nc=1:c
    k1=k1+1;
    j=find(f(:,nc)); 

    nj=numel(j);
    k2=k1+nj-1;
    xb(k1:k2)=nc;
    iif(k1:k2)=f(j,nc);
    if nj==r2  
        j=1:r2;
    else
       dj=[NaN; diff(j)];
       fdj=[1 ;find(dj>1)];
       for nfdj=1:numel(fdj)
           if nfdj~=numel(fdj)
              ij=fdj(nfdj)+(0: find(dj(fdj(nfdj)+1:end)~=1,1,'first')-1); %consecutive
           else
              ij=fdj(nfdj):nj;
           end
           j0=(j(ij(1))* r2/r) ;
           j(ij)=ij-ij(1)+j0;
       end
    end
    yb(k1:k2)=j;
    k1=k2;
end
end
w=1/c;
X=0:w:1-w;
lm=X(xb);

h=1/range(yb);
bm=transf(yb,'R',[0 1-h]);
lm=makeas(iif,1:npf,lm(:));
bm=makeas(iif,1:npf,bm(:));
w=w*.82; h=h*.82;
lm=lm+(w*.18)/2;
else
[y,x]=ind2sub(size(f),find(f));
w=1/c;
h=1/r;
X=0:w:1-w;
Y=0:h:1-h;
lm=X(x); bm=Y(y);
i=f(f>0);
lm=makeas(i,1:npf,lm);
bm=makeas(i,1:npf,bm);
w=w*.82; h=h*.82;

end




% orig=[]; coords=[];
% for nr=1:r
%     j=find(f(nr,:)); 
%     orig=[orig f(j)];
%     j=(j-1)/(c-1)*(1-w);    
%     if numel(j)==c2
%      j=linspace(0,1-w,c2);        
%     else
%     while min(diff(j))<w
%      j=linspace(j(1)*0.95, (1-w-j(1)*0.95)   ,numel(j));
% %      display(num2str(min(diff(j))))
%     end 
%     end
%     coords=[coords; [repmat(nr,[numel(j) 1]) j(:)]];    
% end



% now x is fine, can I gain something in the vertical dimension too?
%let's make it flexible where possible
% fpr=sum(f>1,2); %figures per row
% ad=fpr(1:end-1)+fpr(2:end); %if adjacent rows are  merged
% gain=find(ad<c2); %here they can be merged without exceeding the optimal number of colums
% 
% % try to see if I can insert lower line in upper one...
%     yy=coords(:,1); xx=coords(:,2);
% for nn=1:numel(gain)
%     j=gain(nn);
%     x0=xx(yy==j);
%     x1=xx(yy==j+1);
% end



POS.maximsperfig=npf;
else
    if exist('rt','var')
        %not sure how to do this, it takes to solve x+1/x....
        w=rt*w;
        X=transf(X,'R',1,rt*[w/2 1-w/2]);
        Y=transf(Y,'R',1,[h/2  1-h/2]);     
       lm=X-(w/2);  bm=Y-(h/2); %shifted
    else
        
     w=w*.82; h=h*.82;
        
     
     X=transf(X,'R',1,[w/2 1-w/2]);
     Y=transf(Y,'R',1,[h/2 1-h/2]);   
     

     
    lm=X-(w/2);  bm=Y-(h/2); %shifted
    end
    POS.maximsperfig=numel(lm);
 end
 makeas1=[];
end

w=repmat(w,[npf 1]);
h=repmat(h,[npf 1]);


if exist('space_ax','var') && ~isempty(space_ax) 
hf=(space_ax(3)-space_ax(1));
vf=(space_ax(4)-space_ax(2));
lm=space_ax(1)+(lm*hf);
bm=space_ax(2)+(bm*vf);
w=w*hf;
h=h*vf;
end

if exist('makeas1','var') && ~isempty(makeas1) %only for submode = 1
    ok=true(size(bm));
     for n=1:size(makeas1,1)
         x=makeas1(n,1);
         i=makeas1(n,:);         
         ok(i(2:end))=false;
         bmi=bm(i);  lmi=lm(i);
         bm(x)=min(bmi);
         h(x)=(max(bmi)+h(1) - bm(x));
         lm(x)=min(lmi);
         w(x)=(max(lmi)+w(1) -lm(x));
     end
     w=w(ok);      h=h(ok); 
     lm=lm(ok);     bm=bm(ok);
end

POS.npf=numel(w);     
POS.w=w;
POS.h=h;
POS.bm=bm;
POS.lm=lm;
OUT=POS;   
end

    


