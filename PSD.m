%Fs=1;
psd_data=test_data(400:end,:);
% nfft=2^nextpow2(length(x));
% Pxx=abs(fft(x,nfft)).^2/length(x)/Fs;
% Hpsd=dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',Fs);
% plot(Hpsd)
% set(gca,'xscale','log')

% load parameters.mat;

% load Time.txt;
% t = Time(1:end-1);
% clear Time;



%Intensity=psd_data;
% load Intensity.txt;
%x = sum(Intensity(1:end-1,1:Gturf),2);

% clear Intensity;

%Fs = 1/tsave;
figure
for i=1:1:6
    subplot(2,3,i)
    t=0:0.5:1300;
    x=psd_data(:,i);
    Fs=0.5;
    N = length(t);
    xdft = fft(x);
    xdft = xdft(1:N/2+1);
    psdx = (1/(Fs*N)) * abs(xdft).^2;
    psdx(2:end-1) = 2*psdx(2:end-1);
    freq = 0:Fs/length(x):Fs/2;
    plot(freq,10*log10(psdx), '.');
    grid on
    xlabel('Frequency (Hz)', 'fontsize', 16)
    ylabel('Power/Frequency', 'fontsize', 16)
    set(gca,'xscale','log');
end