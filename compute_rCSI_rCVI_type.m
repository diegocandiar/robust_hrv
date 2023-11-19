function struct_output = compute_rCSI_rCVI_type(ibi, t_ibi, wind, method)

Fs = 4; %4 20
time = t_ibi(1) : 1 / Fs : t_ibi(end);

ibi = ibi(:)';
%% first poincare plot
ibi1 = ibi(1:end-1);
ibi2 = ibi(2:end);

if strcmp(method, 'exact')
    C = cov([ibi1', ibi2']);
    [~, eigenvalues] = eig(C);
    SD01 = sqrt(eigenvalues(1,1));
    SD02 = sqrt(eigenvalues(2,2));
    
    D0 = sqrt( (mean(ibi1))^2 + (mean(ibi2))^2);
end

if strcmp(method, 'approximate')
    sd=diff(ibi); 
    SD01 = sqrt(0.5*std(sd)^2);
    SD02 = sqrt(2*(std(ibi)^2)-(0.5*std(sd)^2));

    D0 = sqrt( (mean(ibi1))^2 + (mean(ibi2))^2);
end

if strcmp(method, 'robust')
    [C, ~] = robustCovHRV([ibi1', ibi2']);
    [~, eigenvalues] = eig(C);
    SD01 = sqrt(eigenvalues(1,1));
    SD02 = sqrt(eigenvalues(2,2));
   
    D0 = sqrt( (trimmean(ibi1,5))^2 + (trimmean(ibi2,5))^2);
end

if strcmp(method, '95%')
    C = robustcov([ibi1', ibi2'],'Method','fmcd','OutlierFraction',0.05);
    [~, eigenvalues] = eig(C);
    SD01 = sqrt(eigenvalues(1,1));
    SD02 = sqrt(eigenvalues(2,2));
   
    D0 = sqrt( (trimmean(ibi1,5))^2 + (trimmean(ibi2,5))^2);
end
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
    t_C(k) = median(t_ibi(ix));

    if strcmp(method, 'exact')
        C = cov([ibi1', ibi2']);
        [~, eigenvalues] = eig(C);
        SD1(k) = sqrt(eigenvalues(1,1));
        SD2(k) = sqrt(eigenvalues(2,2));
        
        D(k) = sqrt( (mean(ibi1))^2 + (mean(ibi2))^2);
    end
    
    if strcmp(method, 'approximate')
        sd=diff(ibi(ix)); 
        SD1(k) = sqrt(0.5*std(sd)^2);
        SD2(k) = sqrt(2*(std(ibi(ix))^2)-(0.5*std(sd)^2));
    
        D(k) = sqrt( (mean(ibi1))^2 + (mean(ibi2))^2);
    end
    
    if strcmp(method, 'robust')
        [C, ~] = robustCovHRV([ibi1', ibi2']);
        [~, eigenvalues] = eig(C);
        SD1(k) = sqrt(eigenvalues(1,1));
        SD2(k) = sqrt(eigenvalues(2,2));
       
        D(k) = sqrt( (trimmean(ibi1,5))^2 + (trimmean(ibi2,5))^2);
    end
    
    if strcmp(method, '95%')
        C = robustcov([ibi1', ibi2'],'Method','fmcd','OutlierFraction',0.05);
        [~, eigenvalues] = eig(C);
        SD1(k) = sqrt(eigenvalues(1,1));
        SD2(k) = sqrt(eigenvalues(2,2));
       
        D(k) = sqrt( (trimmean(ibi1,5))^2 + (trimmean(ibi2,5))^2);
    end


end

SD1 = SD1 - mean(SD1) + SD01;
SD2 = SD2 - mean(SD2) + SD02;
D = D - mean(D) + D0;
% D = 1./D;
% CVI = SD1.*SD2 * 100;
% CSI = SD2./SD1;

CVI = (SD1 * 10 +1); % x 10+1
CSI = SD2 * 1 +1;

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

%% outputs
struct_output = struct;
struct_output.CSI = rCSI;
struct_output.CVI = rCVI; 
struct_output.CSI_HR = Ds_out;
struct_output.CVI_HR = Dv_out; 
struct_output.CSI_HRV = CSI_out;
struct_output.CVI_HRV = CVI_out;
struct_output.time = t_out;
struct_output.fsample = Fs;

end

function fsig = flipsig(sig)

msig = mean(sig);
sig1 =  sig - msig;
sig1 = -sig1;
fsig = sig1 + msig;

end
