f_samp = 330e3;
shoot = 4e3 ;
tolerance = 0.15;
%Band Edge speifications
fs1 = 41e3;
fp1 = 45e3;
fp2 = 65e3;
fs2 = 69e3;

Wc1 = fp1*2*pi/f_samp;
Wc2  = fp2*2*pi/f_samp;
Wer = shoot*2*pi/f_samp;

%Kaiser paramters
A = -20*log10(tolerance);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-8)*165 / (4.57*4*pi));           %empirical formula for N_min

%Window length for Kaiser Window
n=N_min+8 ;
n = 2*n + 1;
%Ideal bandpass impulse response of length "n"
bp_ideal = ideal_lp(Wc2+Wer/2,n) - ideal_lp(Wc1-Wer/2,n);

%Kaiser Window of length "n" help idealwith shape paramter beta calculated above
kaiser_win = (kaiser(n,beta))';

FIR_BandPass = bp_ideal .* kaiser_win;
fvtool(FIR_BandPass,"freq");         %frequency response
legend("magnitude response","phase response");

%magnitude response

figure;
fa = [-(n-1)/2:1:(n-1)/2];
stem(fa, FIR_BandPass);
xlim([-36, 36]);
ylim([-.15, 0.18]);
title("Time domain impulse respone of FIR band pass filter");
xlabel("Samples");

figure;
f = [0:0.01:2.75];
h = freqz(FIR_BandPass,1,f);
plot(330/2/pi*f,abs(h));

ylim([0, 1.2]);
hold on;
set(gca,'XMinorTick','on');
title("Magnitude plot |H(w)| for FIR Band Pass filter");
xlabel('frequency in kHz');
yL = get(gca,'YLim');
line([65 65],yL,'Color','r','Linestyle','--');
xL = get(gca,'XLim');
line(xL,[0.85 0.85],'Color','r','Linestyle','--');
line(xL,[0.15 0.15],'Color','r','Linestyle','--');
line([41 41],yL,'Color','r','Linestyle','--');
line([45 45],yL,'Color','r','Linestyle','--');
line([69 69],yL,'Color','r','Linestyle','--');

