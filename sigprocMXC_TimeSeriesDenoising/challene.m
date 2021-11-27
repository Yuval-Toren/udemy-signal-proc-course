clc;
clear;
load denoising_codeChallenge.mat;
sorted = sort(origSignal);
figure(1)
subplot(4,1,1)
plot(origSignal)
%% median to get rid of all the irregular peaks
med = median(origSignal);

for c = 1:length(origSignal)
    if origSignal(c)<-6 || origSignal(c)>6
        origSignal(c) = med;
    end
end

subplot(4,1,2)
plot(origSignal)

%% run mean smooth
k = 80;
filteredSig = zeros(size(origSignal));

for ii=k+1:length(origSignal)-k-1
    % each point is the average of k surrounding points
    filteredSig(ii) = mean(origSignal(ii-k:ii+k));
end
subplot(4,1,3)
plot(filteredSig(k+1:(length(filteredSig)-k-1)))
subplot(4,1,4)
plot(cleanedSignal)
%% Gauss filter
% fwhm = 50;
% k = 40;
% gtime = -k:k;
% 
% % create Gaussian window
% gauswin = exp( -(4*log(2)*gtime.^2) / fwhm^2 );
% % gauswin = gauswin / sum(gauswin);
% 
% % initialize filtered signal vector
% filtsigG = zeros(size(origSignal));
% 
% % implement the weighted running mean filter
% for i=k+1:length(origSignal)-k-1
%     filtsigG(i) = sum( origSignal(i-k:i+k).*gauswin );
% end
% 
% plot(filtsigG)