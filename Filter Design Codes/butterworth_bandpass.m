%%
%Butterworth Analog LPF parameters
Wc = 1.07;              %cut-off frequency
N = 8;                  %order 

%poles of Butterworth polynomial of degree 8 in the open CLHP 
p1 = Wc*cos(pi/2 + pi/16) + i*Wc*sin(pi/2 + pi/16);
p2 = Wc*cos(pi/2 + pi/16) - i*Wc*sin(pi/2 + pi/16);
p3 = Wc*cos(pi/2 + pi/16+pi/8) + i*Wc*sin(pi/2 + pi/16+pi/8);
p4 = Wc*cos(pi/2 + pi/16+pi/8) - i*Wc*sin(pi/2 + pi/16+pi/8);
p5 = Wc*cos(pi/2 + pi/16+2*pi/8) + i*Wc*sin(pi/2 + pi/16+2*pi/8);
p6 = Wc*cos(pi/2 + pi/16+2*pi/8) - i*Wc*sin(pi/2 + pi/16+2*pi/8);
p7 = Wc*cos(pi/2 + pi/16+3*pi/8) + i*Wc*sin(pi/2 + pi/16+3*pi/8);
p8 = Wc*cos(pi/2 + pi/16+3*pi/8) - i*Wc*sin(pi/2 + pi/16+3*pi/8);

%Band Edge speifications
fp1 = 45;
fs1 = 41;
fs2 = 69;
fp2 = 65;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 330;         
wp1 = tan(fp1/f_samp*pi);
ws1 = tan(fs1/f_samp*pi); 
ws2 = tan(fs2/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);

%Parameters for Bandpass Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8],Wc^N);   %TF with poles p1-p8 and numerator Wc^N and no zeroes
                                                        %numerator chosen to make the DC Gain = 1

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);        %analog LPF Transfer Function
analog_bpf(s) = analog_lpf((s*s + W0*W0)/(B*s));        %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));              %bilinear transformation
discrete_lpf(z) = analog_lpf((z-1)/(z+1));
%coeffs of analog bpf
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bpf
[nz, dz] = numden(discrete_bpf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz,"freq")                                           %frequency response

%magnitude plot (not in log scale) 
f = [0:0.01:2.5];
h = freqz(nz,dz,f);
plot(330/2/pi*f,abs(h));

ylim([0, 1.2]);
hold on;
yL = get(gca,'YLim');
line([41 41],yL,'Color','r','Linestyle','--');
xL = get(gca,'XLim');
line(xL,[0.85 0.85],'Color','r','Linestyle','--');
line(xL,[0.15 0.15],'Color','r','Linestyle','--');
line([45 45],yL,'Color','r','Linestyle','--');
line([65 65],yL,'Color','r','Linestyle','--');
line([69 69],yL,'Color','r','Linestyle','--');
set(gca,'XMinorTick','on');
title("Magnitude plot |H(w)| for Butterworth Band pass filter");
xlabel('Frequency in kHz');

%%