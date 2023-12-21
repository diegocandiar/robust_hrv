function [LF, HF, t_temp] = compute_PWVD(IBI,t_IBI)
Frr = 4; % 4 20
t_temp = t_IBI(1) : 1/Frr : t_IBI(end);
ibi_int = interp1(t_IBI, IBI, t_temp, 'spline'); 
ibi_int = detrend(ibi_int);

Fbins = pow2(12); 
v0 = 0.03; 
tau0 = 0.06; 
lambda = 0.3; 

% Wigner-Ville distribution
ibi_x = hilbert(ibi_int-mean(ibi_int));
ibi_x = ibi_x(:);
n_ibi = 1:length(ibi_x);

[xrow,xcol] = size(ibi_x);
[trow,tcol] = size(n_ibi);

wvx= zeros (Fbins,tcol);  
for icol = 1:tcol
    ti = n_ibi(icol); 
    taumax = min([ti-1,xrow-ti,round(Fbins/2)-1]);
    tau = -taumax:taumax; 
    indices = rem(Fbins+tau,Fbins)+1;
    wvx(indices,icol) = ibi_x(ti+tau,1) .* conj(ibi_x(ti-tau,xcol));
    tau = round(Fbins/2); 
    if ti<=(xrow-tau) && ti>=(tau+1)
        wvx(tau+1,icol) = 0.5 * (ibi_x(ti+tau,1) * conj(ibi_x(ti-tau,xcol))  + ...
                               ibi_x(ti-tau,1) * conj(ibi_x(ti+tau,xcol))) ;
    end
end

wvx= fft(wvx); 
if xcol==1
    wvx=real(wvx); 
end

% Ambiguity function
AF = fft(ifft(wvx).'); clear wvx

if rem(size(AF,2),512)~=0 && rem(size(AF,2),1024)~=0 && rem(size(AF,2),4096)~=0  && rem(size(AF,2),256)~=0
    AF=AF.';
end

[Crow,Ccol] = size(AF);
dy = (-Ccol/2 : Ccol/2-1) / (Ccol/2);

if Crow/2- fix(Crow/2)==0
    dx = (-Crow/2:Crow/2-1) / (Crow/2);
else
    dx = (-Crow/2-1/2:Crow/2-3/2) / (Crow/2);
end

[x,y] = meshgrid(dy,dx);
tau1 = x/tau0;
v1 = y/v0;

beta = 1;
gamma = 1;
alfa = 0;
r = 0;
mu =  (tau1.^2.*(v1.^2).^alfa + (tau1.^2).^alfa .* v1.^2 +2.*r*((tau1.*v1).^beta).^gamma);
clear tau1 v1 x y

% Compute LF HF
t_wind = exp(-pi*( mu.^2).^lambda);
t_wind = t_wind/max(max(t_wind));
Asmooth = (AF.*fftshift(t_wind)).';
TFR_RR = real(ifft(fft(Asmooth).').');

f = linspace(0,Frr/2,Fbins);
[~,LF_ix] = find(f>=0.04 & f<=0.15);
[~,HF_ix] = find(f>=0.15 & f<=0.4);
% 
% [~,LF_ix] = find(f>=0.1 & f<=1);
% [~,HF_ix] = find(f>=1 & f<=3.5);

LF = trapz(TFR_RR(LF_ix,:));
HF = trapz(TFR_RR(HF_ix,:));
end