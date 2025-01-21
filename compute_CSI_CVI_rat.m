function [rCSI, rCVI, Ds_out, Dv_out, CSI_out, CVI_out, t_out] = compute_CSI_CVI_rat(ibi, t_ibi, wind, ks, kp)

Fs = 20; % Hz
time = t_ibi(1) : 1 / Fs : t_ibi(end);

ibi = ibi(:)';
%% first poincare plot
ibi1 = ibi(1:end-1);
ibi2 = ibi(2:end);
C = cov([ibi1', ibi2']);
% C = robustcov([ibi1', ibi2'],'Method','fmcd','OutlierFraction',0.05);
% [C, lam] = robustCovHRV([ibi1', ibi2']);

[eigenvectors, eigenvalues] = eig(C);
SD01 = sqrt(eigenvalues(1,1));
SD02 = sqrt(eigenvalues(2,2));

D0 = sqrt( (trimmean(ibi1,5))^2 + (trimmean(ibi2,5))^2);

%% time varying SD
t1 = time(1);
t2 = t1 +  wind;
ixs = find(t_ibi > t2);
nt = length(ixs)-1;

SD1 = zeros(1,nt);
SD2 = zeros(1,nt);
D = zeros(1,nt);
t_C = zeros(1,nt);

for k = 1 : nt
    i = ixs(k); 

    t2 = t_ibi(i);
    t1 = t_ibi(i)-wind;
    ix = find(t_ibi >= t1 & t_ibi<= t2);
    
    ibi1 = ibi(ix(1:end-1));
    ibi2 = ibi(ix(2:end));
    C = cov([ibi1', ibi2']);
    % C = robustcov([ibi1', ibi2'],'Method','fmcd','OutlierFraction',0.05);
    % [C, lam] = robustCovHRV([ibi1', ibi2']);

    [eigenvectors, eigenvalues] = eig(C);
    SD1(k) = sqrt(eigenvalues(1,1));
    SD2(k) = sqrt(eigenvalues(2,2));
    D(k) = sqrt( (trimmean(ibi1,5))^2 + (trimmean(ibi2,5))^2);
    t_C(k) = median(t_ibi(ix));
end

SD1 = SD1 - mean(SD1) + SD01;
SD2 = SD2 - mean(SD2) + SD02;
D = D - mean(D) + D0;


CVI = SD1 * kp ; 
CSI = SD2 * ks ;

%%
t_out = t_C(1) : 1 / Fs : t_C(end);
CVI_out = interp1(t_C, CVI, t_out, 'Spline');
CSI_out = interp1(t_C, CSI, t_out, 'Spline');
Dv_out = interp1(t_C, D, t_out, 'Spline');
Ds_out = flipsig(Dv_out);

%%

RR_int = interp1(t_ibi, ibi, t_out,'spline');


%%

% rCSI = CSI_out;
% rCVI = CVI_out;

rCSI = Ds_out + CSI_out;
rCVI = Dv_out + CVI_out;

end

function fsig = flipsig(sig)

msig = mean(sig);
sig1 =  sig - msig;
sig1 = -sig1;
fsig = sig1 + msig;

end
