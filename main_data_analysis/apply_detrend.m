function detrended_data = apply_detrend(data, params)
        % helper to apply multiple order detrend on data
        % assume data structure as before: each column a subject, evolution
        % of time series along first dimension (rows)

        detrended_data = nan(size(data));
        n_subj = size(data,2);
        
        
        for iPoly = 1:n_subj

            vect_data = data(:, iPoly);
            current_poly = polyfit(params.time_bins, vect_data, params.detrend_flag);
            trend = polyval(current_poly, params.time_bins);
            detrended_data(:, iPoly) = vect_data-trend;

        end

    end