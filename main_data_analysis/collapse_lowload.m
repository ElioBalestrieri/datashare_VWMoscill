

load('../RAW_data_collapsed.mat')
n_subject = size(mat_data, 3);
    
HR_L02 = nan(16, n_subject);

want_all_trials = true;



for iSubj = 1:n_subject

    curr_mat = mat_data(:,:,iSubj);
        
    % warning --nan + get delta values 
    delta_vals = unique(curr_mat(4,:));
    
    memAcc_Lmask = curr_mat(2,:)==curr_mat(3,:);
    memAcc_Lmask(curr_mat(1,:)==0)=1; % set to correct trials with 0 load
    
    if want_all_trials 
        
        memAcc_Lmask = true(size(memAcc_Lmask)); % take all trials into acount regardless of memory outcome
        
    end
    
    loadmask = (curr_mat(1, :) == 0) | (curr_mat(1, :) == 2);
    

    loopDelta = 0;
    for iDelta = delta_vals

        loopDelta = loopDelta+1;

        % find current delta t in data
        deltaT_Lmask = curr_mat(4,:)==iDelta;

        % mask for current load, delta t in right memory trials
        delta_load_mem_Lmask = deltaT_Lmask & memAcc_Lmask & loadmask;

        % mask for HR (only flash present)
        flash_present = curr_mat(6,:)==1;
        HR_lMask = delta_load_mem_Lmask&flash_present;

        % put HR into matrix, for each load, each delta, each
        % participant
        curr_HR = nanmean(curr_mat(7,HR_lMask));
        HR_L02(loopDelta, iSubj) = curr_HR;

        
    end
    
end

if want_all_trials

    save('../HR_allsubj_L02_alltrls.mat','HR_L02')

else

    save('../HR_allsubj_L02.mat','HR_L02')

end
