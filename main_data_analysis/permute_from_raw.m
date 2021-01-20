function HR_mat_perm = permute_from_raw(want_all_trls, n_iter)

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

    HR_mat_perm = nan(3, 16, n_subj, n_iter);

    for iSubj = 1:n_subj

        curr_mat = mat_data(:,:,iSubj);


        if want_all_trls

            vect_acc_mem = true(1,1440);

        else

            % not take into consideration wrong responses
            vect_acc_mem = curr_mat(2,:)==curr_mat(3,:);

        end


        loopLoad = 0;
        for iLoad = [0 2 4]

            loopLoad = loopLoad+1;
            vect_curr_load = curr_mat(1,:)==iLoad;

            % correct -every mem trial in load0 is right
            if iLoad == 0
                vect_acc_mem(vect_curr_load) = 1;
            end

            % logical mask gathering current load and memory
            lgcl_mask_load_mem_acc = vect_curr_load&vect_acc_mem;

            % create reduced matrix containing only flash_presence, flash
            % response and delta t        
            reduced_mat = curr_mat([4 6 7], lgcl_mask_load_mem_acc);
            vect_deltaT = reduced_mat(1,:); length_vect_DT = numel(vect_deltaT); 

            % start n iteration of datasmpling

            swap_reduced_mat = reduced_mat;
            for iIter = 1:n_iter
                
                vect_deltaT = vect_deltaT(randsample(1:length_vect_DT,...
                    length_vect_DT));

                swap_reduced_mat(1,:,:) = vect_deltaT(randsample(...
                    1:length_vect_DT, length_vect_DT));

                loopDelta = 0;
                for iDelta = delta_ts

                    loopDelta = loopDelta +1;
                    mini_mat = swap_reduced_mat([2 3], swap_reduced_mat(1,:)...
                        ==iDelta);

                    vect_HR = mini_mat(1,:)&mini_mat(2,:);

                    HR = sum(vect_HR)/sum(mini_mat(1,:));

                    HR_mat_perm(loopLoad,loopDelta,iSubj, iIter) = HR;

                end

            end       

        end

    end

end