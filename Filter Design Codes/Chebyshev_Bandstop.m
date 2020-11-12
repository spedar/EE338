%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 4;

% Open CLHP Poles of the Chebyshev Polynomial of order 4
p1 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p2 = -sin(pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(pi/(2*N))*cosh(asinh(1/epsilon)/N);
p3 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)+i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);
p4 = -sin(3*pi/(2*N))*sinh(asinh(1/epsilon)/N)-i*cos(3*pi/(2*N))*cosh(asinh(1/epsilon)/N);        

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2);          %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+epsilon*epsilon))];        % even order, DC Gain set as 1/(1+ epsilon^2)^0.5

%Band Edge speifications
fs1 = 39.2;
fp1 = 35.2;
fp2 = 63.2;
fs2 = 59.2;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 260;
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));     %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));          %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz,"freq");                                       %frequency response in dB

%magnitude plot (not in log scale) 
f = [0:0.01:2.5];
h = freqz(nz,dz,f);
plot(260/2/pi*f,abs(h));

ylim([0, 1.2]);
hold on;
yL = get(gca,'YLim');
line([35 35],yL,'Color','r','Linestyle','--');
xL = get(gca,'XLim');
line(xL,[0.85 0.85],'Color','r','Linestyle','--');
line(xL,[0.15 0.15],'Color','r','Linestyle','--');
line([39 39],yL,'Color','r','Linestyle','--');
line([59.2 59.2],yL,'Color','r','Linestyle','--');
line([63.2 63.2],yL,'Color','r','Linestyle','--');
set(gca,'XMinorTick','on');
title("Magnitude plot |H(w)| for Chebyshev Band stop filter");
xlabel('frequency in kHz');
grid