function [P,f,Y]=PSD_FFT(y,sr,dim,psdfft_sw,L,winkind,tukey_opt)
%        [P,f,Y]=PSD_FFT(y,sr,dim,psdfft_sw,L,winkind,tukey_opt)
%psdfft_sw 1 for PSD 
%          2 for fft
% PSD (Power Spectral Density, the power carried by the signal per unit
% frequency) is computed with the Welch's method. Default Welch's
% segmentation (i.e. 8 sections with 50% overlap) is overriden by L
%and divides y in sections L long, individually windowed with winkind and
%with 50% overlapping. Note that if L equals size(y,dim) no segmentation is
%performed: if a window is used ('winkind' string) it behaves as a
%periodogram, otherwise it is just a raw squared fft amplitude estimate.
% General remarks:
% -PSD for voltage signals it is customary to use units of V2/Hz. Power can
% be the actual physical power, or for convenience with
% abstract signals, can be defined as the squared value of the signal.
% Multiplied by time it gives Energy in the generalized sense of signal
% processing i.e. variance of the signal (V2s2). Divided by impedance
% it gives energy in the physical sense (V2S2/Ohm, that is J/Hz). 
%-computing a periodogram from a finite-lengthdigital sequence using the
%fast Fourier transform (FFT) is not a good spectral estimate because of:
%       SPECTRAL BIAS (the spectral bias problem arises from a sharp
%truncation of the sequence <-- multiply the finite sequence by a window function
%       each frequency VARIANCE DOES NOT DECREASE as the
% number of samples used in the computation increases <--smoothing the
% periodogram
%- Welch's method is an approach to spectral density estimation. It is an
% improvement on the standard periodogram spectrum estimating method and on
% Bartlett's method. It reduces noise in the estimated power
% spectra  in exchange for reducing the frequency resolution. Due to the
% noise caused by imperfect and finite data, the noise reduction from
% Welch's method is often desired.
%Choosing different windows allows tradeoffs between resolution (e.g.,
%using a rectangular window) and sidelobe attenuation (e.g., using a Hann
%window). An Hamming window has a 42.5 dB sidelobe attenuation. This may
%mask spectral content below this value (relative to the peak spectral
%content)....
% TUKEY windows are cosine-tapered windows:
% 'tukey_opt' is the ratio of taper to constant sections
%             it ranges from 0 (=rectwin) to 1  (=hann)
%
%      [P,f,Y]=PSD_FFT(y,sr,dim,psdfft_sw,L,winkind,tukey_opt)

% "Just as the Power Spectral Density (PSD) is the Fourier transform of the
% auto-covariance function we may define the Cross Spectral Density (CSD)
% as the Fourier transform of the cross-covariance function."
%        [P,f,Y]=PSD_FFT(y,sr,dim,psdfft_sw,L,winkind)







if ~nargin
    help PSD_FFT
else
    
if ~exist('dim','var')||isempty(dim)
    [~,dim]=max(size(y));
end
if ~exist('psdfft_sw','var')||isempty(psdfft_sw)
    psdfft_sw=1;
end








if psdfft_sw~=2 %periodogram, Welch's method

if exist('L','var') && ~isempty(L)
if  exist('winkind','var') && ~isempty(winkind) && ~exist('tukey_opt','var')
    wn=['window(@' winkind ',L)'];
    try
        eval(wn);
    catch
        display('Window type is not recognized, defualt is used [hamming]')
        wn='window(@hamming,L)';
    end
elseif exist('tukey_opt','var')
    wn='tukeywin(L,tukey_opt)';
else
    wn='window(@hamming,L)';
end

else %no window if length is not provided
    wn='[]';
end
wn=eval(wn);

if dim==2
[~,f] =  pwelch(y(1,:,1),numel(wn),[],[],sr);
P=zeros(size(y,1),numel(f),size(y,3));
for d3=1:size(y,3)
for d2=1:size(y,1)
 P(d2,:,d3) =  pwelch(y(d2,:,d3),wn,[],[],sr);
end
end
else
[~,f] =  pwelch(y(:,1,1),numel(wn),[],[],sr);
P=zeros(numel(f),size(y,2),size(y,3));
for d3=1:size(y,3)
for d2=1:size(y,2)
 P(:,d2,d3) = pwelch(y(:,d2,d3),wn,[],[],sr);
end
end     
end

 
else %FFT
Y=fft(y,[],dim);
L=size(Y,dim);

if exist('sr','var')||~isempty(sr)
    f = (sr/2*linspace(0,1,L/2+1))';
%     P=abs(Y).^2; 
    P=abs(Y).^2/L/sr; %normalized

    s2=size(P);     s=ones(3,1);    s(1:numel(s2))=s2;  s(dim)=numel(f);
    P=P(1:s(1),1:s(2),1:s(3));
    Y=abs(Y(1:s(1),1:s(2),1:s(3)));
end 
end
end










