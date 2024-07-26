function [E_out, Arr_out, Acc_out, t_out] = compute_SOPP(ibi, t_ibi, wind)

Fs = 4;
time = t_ibi(1) : 1 / Fs : t_ibi(end);

%% time varying SD
t1 = time(1);
t2 = t1 +  wind;
ixs = find(t_ibi > t2);
nt = length(ixs)-1;

E = zeros(1,nt);
Arr = zeros(1,nt);
Acc = zeros(1,nt);
t_C = zeros(1,nt);

for k = 1 : nt
    i = ixs(k); 

    t2 = t_ibi(i);
    t1 = t_ibi(i)-wind;
    ix = find(t_ibi >= t1 & t_ibi<= t2);
    
    ibi1 = ibi(ix(2:end-1)) - ibi(ix(1:end-2));
    ibi2 = ibi(ix(3:end)) - ibi(ix(2:end-1)); 

    % Initialize lists to hold points in each quadrant
    Q1 = [];
    Q2 = [];
    Q3 = [];
    Q4 = [];
    
    % Categorize points into quadrants
    for i = 1:length(ibi1)
        x = ibi1(i);
        y = ibi2(i);
        
        if x >= 0 && y >= 0
            Q1 = [Q1; x, y];
        elseif x < 0 && y > 0
            Q2 = [Q2; x, y];
        elseif x <= 0 && y <= 0
            Q3 = [Q3; x, y];
        elseif x > 0 && y < 0
            Q4 = [Q4; x, y];
        end
    end

    Q1_var = calculate_var(Q1);
    Q2_var = calculate_var(Q2);
    Q3_var = calculate_var(Q3);
    Q4_var = calculate_var(Q4);

    Arr(k) = Q2_var-Q4_var; % 
    Acc(k) = Q3_var-Q1_var; %  

    E0 = 0;
    Q = [length(Q1) length(Q2) length(Q3) length(Q4)];
    proportions = Q / sum(Q);
    % Compute entropy
    for i = 1:length(proportions)
        if proportions(i) > 0
            E0 = E0 - proportions(i) * log2(proportions(i));
        end
    end

    E(k) = E0;

    t_C(k) = median(t_ibi(ix));
end


%%
t_out = t_C(1) : 1 / Fs : t_C(end);
E_out = interp1(t_C, E, t_out, 'Spline');
Arr_out = interp1(t_C, Arr, t_out, 'Spline');
Acc_out = interp1(t_C, Acc, t_out, 'Spline');

end


function varq = calculate_var(quadrant_points)
    if isempty(quadrant_points)
        varq = 0; % Handle empty quadrants
    else
        distances = sqrt(quadrant_points(:,1).^2 + quadrant_points(:,2).^2);
        varq = 2 * std(distances);
    end
end

