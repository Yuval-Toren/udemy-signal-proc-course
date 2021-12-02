clc;clear;
load spectral_codeChallenge.mat;
n = length(signal);
spow = abs(fft(signal)/n).^2;
hz = linspace(0,srate/2,floor(n/2)+1);
plot(hz,spow(1:length(hz)))
%% ref
% [powspect,frex,time] = spectrogram(signal,hann(1000),100,[],srate);
% 
% % Octave uses the following line instead of the above line
% %[powspect,frex,time] = specgram(detrend(bc(:,2)),1000,fs,hann(1000));
% 
% 
% % show the time-frequency power plot
% figure(1), clf
% imagesc(time,frex,abs(powspect).^2)
% axis xy
% xlabel('Time (sec.)'), ylabel('Frequency (Hz)')
% set(gca,'clim',[0 1]*2,'ylim',frex([1 dsearchn(frex(:),15000)]),'xlim',time([1 end]))
% colormap hot
%% code
% window length in seconds*srate
winlength = (1)*srate;

% number of points of overlap
nOverlap = round(srate/2);

% window onset times
winonsets = 1:winlength-nOverlap:n-winlength;

% note: different-length signal needs a different-length Hz vector
hzW = linspace(0,srate/2,floor(winlength/2)+1);

% Hann window
hannw = .5 - cos(2*pi*linspace(0,1,winlength))./2;

% initialize the power matrix (windows x frequencies)
% signalW = zeros(1,length(hzW));
signalW = [];

% loop over frequencies
for wi=1:length(winonsets)
    
    % get a chunk of data from this time window
    datachunk = signal(winonsets(wi):winonsets(wi)+winlength-1);
    
    % apply Hann taper to data
    datachunk = datachunk .* hannw;
    
    % compute its power
    tmppow = abs(fft(datachunk)/winlength).^2;
    
    % enter into matrix
    signalW = [signalW, tmppow'];
end

% divide by N
% signalW = signalW / length(winonsets);
signalW = signalW(1:length(hzW),:);
imagesc(time/srate,hz,signalW)
axis xy
xlabel('Time (sec.)'), ylabel('Frequency (Hz)')
colormap hot