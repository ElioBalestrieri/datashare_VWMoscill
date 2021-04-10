function [detrended_data, trends_all_subj, coeff_polynomials] = ...
    do_detrend(data, n_subj, params)
    % helper to apply multiple order detrend on data
    % assume data structure as before: each column a subject, evolution
    % of time series along first dimension (rows)

    [trends_all_subj, detrended_data] = deal(nan(size(data)));
    
    coeff_polynomials = nan(n_subj, params.detrend_flag+1);
    

    for iPoly = 1:n_subj

        vect_data = data(:, iPoly);
        current_poly = polyfit(params.time_bins, vect_data, params.detrend_flag);
        coeff_polynomials(iPoly,:) = current_poly;
        trend = polyval(current_poly, params.time_bins);
        trends_all_subj(:, iPoly) = trend;
        detrended_data(:, iPoly) = vect_data-trend;
 
    end

%     vect_data = data;
%     current_poly = polyfit(params.time_bins, vect_data, params.detrend_flag);
%     trend = polyval(current_poly, params.time_bins);
%     detrended_data = vect_data-trend;

    
end
