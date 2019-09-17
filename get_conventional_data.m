function [frequency_inner frequency_outer congruent correct] = get_conventional_data(run_idx, present_inner, present_outer, Results)
for trial = 1:length(Results.PDir)
    frequency_inner(trial,1) = length(find(diff(Results.discrete_steps{trial})...
        ~= 0)) / length(Results.discrete_steps{trial});
    
    if Results.present_outer == 0
        frequency_inner(trial,1) = length(find(diff(Results.discrete_steps{trial})...
            ~= 0)) / length(Results.discrete_steps{trial});
        frequency_outer(trial,1) = NaN;
        congruent(trial, 1) = NaN;
        correct(trial, 1) = NaN;

    else
        frequency_outer(trial,1) =  length(find(diff(Results.template.discrete_steps{trial})...
            ~= 0)) / length(Results.template.discrete_steps{trial});
        
        if present_inner == 0
            frequency_inner(trial,1) = NaN;
            congruent(trial, 1) = NaN;
            
            if length(Results.discrete_steps{trial}) == 1
            correct(trial, 1) = Results.random_vector(trial) + Results.discrete_steps{trial};
            else
               correct(trial, 1) = NaN;  
                
            end
 
        else
            frequency_inner(trial,1) = length(find(diff(Results.discrete_steps{trial}...
                .* Results.template.discrete_steps{trial} .*2)~=0)) /...
                length(Results.discrete_steps{trial}); 
            correct(trial, 1) = NaN;
            congruent(trial, 1) = Results.discrete_steps{trial}; 
        end
    end 
end

congruent(congruent == 2) = 1;
congruent(congruent == -2) = 0;

correct(find(abs(correct) == 2)) = 1;

end