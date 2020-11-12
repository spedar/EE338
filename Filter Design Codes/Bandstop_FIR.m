f_samp = 260e3;
shoot = 4e3 ;
tolerance = 0.15;

%Band Edge speifications
fs1 = 39.2e3;
fp1 = 35.2e3;
fp2 = 63.2e3;
fs2 = 59.2e3;

Wer = shoot*2*pi/f_samp;
Wc1 = fs1*2*pi/f_samp;
Wc2  = fs2*2*pi/f_samp;
%Kaiser paramters
A = -20*log10(tolerance);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

Wn = [(fs1+fp1)/2 (fs2+fp2)/2]*2/f_samp;        %average value of the two paramters
N_min = ceil((A-8) / (4.57*Wer));       %empirical formula for N_min

%Window length for Kaiser Window
n=N_min  ;
n= 2*n+1 ;

%Ideal bandstop impulse response of length "n"

bs_ideal =  ideal_lp(pi,n) -ideal_lp(Wc2+Wer/2,n) + ideal_lp(Wc1-Wer/2,n);

%Kaiser Window of length "n" with shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandStop = bs_ideal .* kaiser_win;
fvtool(FIR_BandStop,"freq");         %frequency response
legend("magnitude response","phase response");

%magnitude response
figure;
fa = [-(n-1)/2:1:(n-1)/2];
stem(fa, FIR_BandStop);
xlim([-23, 23]);
ylim([-.15, 0.83]);
title("Time domain impulse respone of FIR band stop filter");
xlabel("Samples");

figure;
f = [0:0.01:2.75];
h = freqz(FIR_BandStop,1,f);
plot(f_samp/1e3/2/pi*f,abs(h));

ylim([0, 1.2]);
hold on;
set(gca,'XMinorTick','on');
title("Magnitude plot |H(w)| for FIR Band Stop filter");
xlabel('frequency in kHz');
yL = get(gca,'YLim');
line([39.2 39.2],yL,'Color','r','Linestyle','--');
xL = get(gca,'XLim');
line(xL,[0.85 0.85],'Color','r','Linestyle','--');
line(xL,[0.15 0.15],'Color','r','Linestyle','--');
line([59.2 59.2],yL,'Color','r','Linestyle','--');
line([35.2 35.2],yL,'Color','r','Linestyle','--');
line([63.2 63.2],yL,'Color','r','Linestyle','--');
