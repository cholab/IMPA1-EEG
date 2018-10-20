function [M,i]=stitch(M,stitchmode, detrending, p, overlap, avrange, fq,sr)
%To adjust a 3D (el x tp x trials) dataset M to optimize stitching (i.e. improve frequency analysis
%at boundaries):
%            [Ms,i]=stitch(M,stitchmode, detrending)
%Optionally, detrending is applied beforehand (type stitch(M,0,1) to run
%detrending only).
%Stitching is implemented through 4 different strategies (see stitchmodes, [2]).
%To enter optional inputs for padding:
%            [Ms,i]=stitch(M,stitchmode, detrending, p, overlap, avrange, fq,sr)
%        p         is pad length [100]
%        For EXPAD:
%        overlap   is the overlap between extrapolated and real data  for smooth transition,
%                  note that it defines also the length of  boundaries for which
%                  phase/amplitude information is discarded [10]
%        avrange   small window where average amplitude is calculated for initial extrapolation
%        fq        frequencies used for phase-amplitude extrapolation fq=[1:5 6 8 10 20 30 40]; 
%        sr        sampling rate [250]
%STITCH MODES are:
%  1 zero step - no padding, segment n offset (segment(n,:) - (segment(n,1)- segment(n-1,end)) ) 
%  2 padding: linear - a pad long 'p' is added that linearly joins segment(n-1,end)-segment(n,1)
%  3 padding: anticipate step - the pad move the (n-1)-(n) step at a distance p/2
%  4 padding: extrapolation based smart padding (EX.PAD) 
%             -phase PH for fq is extrapolated at boundaries (discarding 'overlap' points)
%             -synthetic data PF*A are created, amplitude A is 0 at the
%             extremities and linearly rise to the mean value of 'avrange'
%             real data 
%             -smooth transition from real data to synthetic is implemented
%The vector of real data is 'i'.


if ~nargin
    help stitch
    return
end

if ~exist('stitchmode','var')||isempty(stitchmode)
    stitchmode=2;
end
if ~exist('detrending','var')||isempty(detrending)
    detrending=1;
end
if ~exist('p','var')||isempty(p)
    p=100;
end
if ~exist('overlap','var')||isempty(overlap)
    overlap=10;
end
if ~exist('avrange','var')||isempty(avrange)
    avrange=40;
end
if ~exist('fq','var')||isempty(fq)
    fq=[1:5 6 8 10 20 30 40]; 
end
if ~exist('sr','var')||isempty(sr)
    sr=250; 
end


[d1 d2 d3]=size(M);
if detrending
    fprintf('Detrending...')
    for n=1:d3
           M(:,:,n)=detrend(M(:,:,n)')';            
    end  
end

switch stitchmode
    case 1 % zero step, no padding
        fprintf('Zero step stitching, no padding...')
        for n=2:d3
           steps=(M(:,end,n-1) - M(:,1,n));
           M(:,:,n)=M(:,:,n)+repmat(steps,[1 d2]); 
        end    
        i=true(1,d2);     
    case 2 %linear padding
        fprintf('Stitching with linear padding...')        
        Ms=zeros(d1,p+d2,d3);
        ir=linspace(0,1,p);
        i=[false(1,p) true(1,d2)];
        Ms(:,:,1)=[repmat(M(:,1,1),[1 p]) M(:,:,1) ];                 
        for n=2:d3
           pad=[M(:,end,n-1)  M(:,1,n)];
           pad=interp1([0 1], pad', ir);
           Ms(:,:,n)=[pad' M(:,:,n) ]; 
        end   
        M=Ms;
    case 3 %padding with anticipated step
        fprintf('Stitching with padding and anticipated step...')
        Ms=zeros(d1,p+d2,d3);
        ir1=round(p/2); ir2=p-ir1;
        Ms(:,:,1)=[repmat(M(:,1,1),[1 p]) M(:,:,1) ];          
        for n=2:d3
           pad=[repmat(M(:,end,n-1),[1 ir1]) repmat(M(:,1,n),[1 ir2])];
           Ms(:,:,n)=[pad M(:,:,n) ]; 
        end           
        i=[false(1,p) true(1,d2)];
        M=Ms;        
    case 4
        %EXPAD
        fprintf('EX.PAD stitching...')
        p=p*2;
        Ms=zeros(d1,p*2+d2,d3);
        nf=numel(fq);
        hl=round(d2/2+p);
        w=sin([ones(1,p)*pi/2  interp1([0 1],[pi/2 -pi/2] ,linspace(0,1,overlap))   zeros(1,d2-overlap*2)-pi/2   interp1([0 1],[-pi/2 pi/2] ,linspace(0,1,overlap))  ones(1,p)*pi/2]);
        w=(w/2)+.5;
        w=[w;1-w];
        for n=1:d3
            display(['trial ' num2str(n) '/' num2str(d3)])
            x=[zeros(d1,p) M(:,:,n) zeros(d1,p)];
            W=wavtransform(M(:,:,n),fq,sr,0,2);
            for m=1:d1
               PH=interp1(overlap+1:d2-overlap,angle(W(:,overlap+1:d2-overlap,m))',1-p:d2+p,'linear','extrap')';
               PH=wrapToPi(PH);
               A=abs(W(:,overlap+1:d2-overlap,m));
%                A1=interp1([0 1],[zeros(nf,1) mean(A(:,1:avrange),2)]', linspace(0,1,p+overlap))';
%                A2=interp1([0 1],[ mean(A(:,d2-2*overlap-avrange+1:end),2) zeros(nf,1)]', linspace(0,1,p+overlap))';                
               A1=repmat(mean(A(:,1:avrange),2),[1 p+overlap]);
               A2=repmat(mean(A(:,d2-2*overlap-avrange+1:end),2),[1 p+overlap]);
       
               synth=sum(cos(PH).*[A1 A A2]);
               fc=range(x(m,:))/range(synth);
               synth=synth*fc;               
               synth1=synth(1:hl)-(synth(p+overlap)-x(m,p+overlap));
               synth2=synth(hl+1:end)-(synth(end-p-overlap)-x(m,end-p-overlap));
               Ms(m,:,n)=sum((([[synth1 synth2]; x(m,:)].*w)));               
            end 
        end  
        w=sin([interp1([0 1],[pi/2 -pi/2] ,linspace(0,1,p))]);        
        w=(w/2)+.5;
        p=p/2;
        w=repmat(w(1:p),[d1 1]);
        w(:,:,2)=1-w;
        Ms2=Ms(:,p+1:end-p,:);
        r=d2+3*p+1;
        for n=1:d3
            if n>1
            A=Ms(:,:,n-1);
            B=Ms(:,:,n);
            pad(:,:,1)=A(:,r:end);
            pad(:,:,2)=B(:,1:p);
            pad=sum(pad.*w,3);
            Ms2(:,1:p,n)=pad;
            end
            if n<d3
            A=Ms(:,:,n);
            B=Ms(:,:,n+1);
            pad(:,:,1)=A(:,r:end);
            pad(:,:,2)=B(:,1:p);
            pad=sum(pad.*w,3);
            Ms2(:,end-p+1:end,n)=pad;
            end
        end
        i=round(overlap/2);      
        i=[false(1,p+i) true(1,d2-i*2) false(1,p+i) ];
        M=Ms2;    
end

fprintf(' done.\n')
