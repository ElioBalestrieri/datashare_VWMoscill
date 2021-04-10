function VWM_permuted = VWM_permute_from_raw(n_iter)

    load('../RAW_data_collapsed.mat');

    %% reminder -data structure-
    % 1 row: load condition
    % 2 row: mem equality condition
    % 3 row: subj mem response
    % 4 row: deltaT
    % 5 row: timestamp delta
    % 6 row: flash_presence
    % 7 row: subj flash response

    delta_ts = unique(mat_data(4,:,:))';    
    
    n_subj = size(mat_data,3);

    VWM_permuted = nan(2, 16, n_subj, n_iter);

    for iSubj = 1:n_subj

        curr_mat = mat_data(:,:,iSubj);

        loopLoad = 0;
        for iLoad = [2 4]

            loopLoad = loopLoad+1;
            vect_curr_load = curr_mat(1,:)==iLoad;
            
            % count only flash present trials
            isflashpresent = curr_mat(6, :)==1;

            % logical mask gathering current load and memory
            lgcl_mask_load_flashon = vect_curr_load & isflashpresent;

            % create reduced matrix containing only flash_presence, flash
            % response and delta t        
            reduced_mat = curr_mat([2 3 4], lgcl_mask_load_flashon);
            length_vect_DT = size(reduced_mat, 2); 

            % start n iteration of datasampling
            
            for iIter = 1:n_iter
                
                swap_reduced_mat = reduced_mat;
                
                vect_deltaT = reduced_mat(3,:);
                
                vect_deltaT = vect_deltaT(randsample(1:length_vect_DT,...
                    length_vect_DT));

                swap_reduced_mat(3,:) = vect_deltaT;

                loopDelta = 0;
                for iDelta = delta_ts

                    loopDelta = loopDelta +1;
                    mini_mat = swap_reduced_mat([1, 2], swap_reduced_mat(3,:)...
                        ==iDelta);
                
                    vect_acc = mini_mat(2, :) == mini_mat(1, :);
                    
                    VWM_permuted(loopLoad,loopDelta,iSubj, iIter) = mean(vect_acc);

                end

            end       

        end

    end

end