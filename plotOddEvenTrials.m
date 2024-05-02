  % Sara Molas Medina
% 20th February 2023

% Plot raw spiked and positional

function [ratemaps, ratemap_odd, ratemap_even] = plotOddEvenTrials (ratemaps, in, onlyPlaceCells, plot_figure)

counterPlaceCells = [0, 0, 0];
 cellIdx = NaN (size(ratemaps.cellN,1),3);
    
    for s = 1:2:3
        if in.whichSessionsToUse(s) == 1
            for c = 1:length(ratemaps.cellN)
                if onlyPlaceCells == 1
                    if in.whichSessionsToUse(1)== 1 && ratemaps.cellN{c}.session(1).Track.cellOn == 0 && in.whichSessionsToUse(3)== 1 && ratemaps.cellN{c}.session(3).Track.cellOn == 0 
                        continue,
                    elseif in.whichSessionsToUse(1)== 1 && ratemaps.cellN{c}.session(1).Track.cellOn == 0 && in.whichSessionsToUse(3)== 0
                        continue,
                    elseif in.whichSessionsToUse(1)== 0  && in.whichSessionsToUse(3)== 1 && ratemaps.cellN{c}.session(3).Track.cellOn == 0
                        continue, 
                    else 
                        counterPlaceCells(s) = counterPlaceCells(s) +1;
                        cellIdx(c,s) = c;
                         vector_binSpks_even = zeros(1, size (ratemaps.cellN{c,1}.session(s).Track.trials{1,1}.ratemap,2));
                         vector_binSpks_odd = zeros(1, size (ratemaps.cellN{c,1}.session(s).Track.trials{1,1}.ratemap,2));
                         vector_binPos_even = zeros(1, size (ratemaps.cellN{c,1}.session(s).Track.trials{1,1}.ratemap,2));
                         vector_binPos_odd = zeros(1, size (ratemaps.cellN{c,1}.session(s).Track.trials{1,1}.ratemap,2));

                        for trial_i = 1:length (ratemaps.cellN{c,1}.session(s).Track.trials)

                            if mod(trial_i,2) == 0 %is even trial
                                binned_spikes = ratemaps.cellN{c,1}.session(s).Track.trials{1,trial_i}.binSpks;
                                binned_pos = ratemaps.vr(s).Track.trial{trial_i}.posBinData;
                                vector_binSpks_even = binned_spikes + vector_binSpks_even;
                                vector_binPos_even =  binned_pos + vector_binPos_even;
                            else % is odd trial
                                binned_spikes = ratemaps.cellN{c,1}.session(s).Track.trials{1,trial_i}.binSpks;
                                binned_pos = ratemaps.vr(s).Track.trial{trial_i}.posBinData;
                                vector_binSpks_odd = binned_spikes + vector_binSpks_odd;
                                vector_binPos_odd =  binned_pos + vector_binPos_odd;

                            end 
                        end 
                        %fprintf('cell %d session %d cell_id %d \n',c,s,ratemaps.cellN{c,1}.session(s).cellInd(1,2));
                        ratemap_odd.session(s).unsmoothed_ratemap(counterPlaceCells(s),:)= vector_binSpks_odd ./ vector_binPos_odd;
                        ratemap_even.session(s).unsmoothed_ratemap(counterPlaceCells(s),:)= vector_binSpks_even ./ vector_binPos_even;
                        visited_pos_odd = ones(1, size(vector_binPos_odd, 2));
                        unvisited_pos_odd = find(vector_binPos_odd==0);
                        visited_pos_odd (1, unvisited_pos_odd) = 0; 
                        ratemap_odd.session(s).visitedPos(counterPlaceCells(s),:) = logical(visited_pos_odd); 

                        visited_pos_even = ones(1, size(vector_binPos_even, 2));
                        unvisited_pos_even = find(vector_binPos_even==0);
                        visited_pos_even (1, unvisited_pos_even) = 0; 
                        ratemap_even.session(s).visitedPos(counterPlaceCells(s),:) = logical(visited_pos_even); 

                        ratemap_odd.session(s).ratemap(counterPlaceCells(s),:) = make_smooth_ratemap(ratemap_odd.session(s).visitedPos(counterPlaceCells(s),:), ratemap_odd.session(s).unsmoothed_ratemap(counterPlaceCells(s),:), in.rmSmooth, 'gaus', 'norm'); 
                        ratemap_even.session(s).ratemap(counterPlaceCells(s),:) = make_smooth_ratemap(ratemap_even.session(s).visitedPos(counterPlaceCells(s),:), ratemap_even.session(s).unsmoothed_ratemap(counterPlaceCells(s),:), in.rmSmooth, 'gaus', 'norm');  
                        ratemap_odd.session(s).peakR(counterPlaceCells(s),:)= max(ratemap_odd.session(s).ratemap(counterPlaceCells(s),:));
                        ratemap_even.session(s).peakR(counterPlaceCells(s),:)= max(ratemap_even.session(s).ratemap(counterPlaceCells(s),:));
                        ratemap_odd.session(s).ratemap_notVisitedVector(counterPlaceCells(s),:) = make_smooth_ratemap(vector_binPos_odd, vector_binSpks_odd, in.rmSmooth, 'gaus', 'norm'); 
                        ratemap_even.session(s).ratemap_notVisitedVector(counterPlaceCells(s),:) = make_smooth_ratemap(vector_binPos_even, vector_binSpks_even, in.rmSmooth, 'gaus', 'norm');  
                        ratemap_odd.session(s).norm_ratemap(counterPlaceCells(s),:) = ratemap_odd.session(s).ratemap(counterPlaceCells(s),:) / max(ratemap_odd.session(s).ratemap(counterPlaceCells(s),:));
                        ratemap_even.session(s).norm_ratemap (counterPlaceCells(s),:) = ratemap_even.session(s).ratemap(counterPlaceCells(s),:) / max(ratemap_even.session(s).ratemap(counterPlaceCells(s),:));
                        clear vector_binSpks_even vector_binSpks_odd 
                        
                        ratemaps.cellN{c,1}.session(s).norm_ratemap_odd = ratemap_odd.session(s).norm_ratemap(counterPlaceCells(s),:);
                        ratemaps.cellN{c,1}.session(s).ratemap_odd = ratemap_odd.session(s).ratemap(counterPlaceCells(s),:);
                        ratemaps.cellN{c,1}.session(s).peakR_odd = ratemap_odd.session(s).peakR(counterPlaceCells(s),:);
                        ratemaps.cellN{c,1}.session(s).norm_ratemap_even = ratemap_even.session(s).norm_ratemap(counterPlaceCells(s),:);
                        ratemaps.cellN{c,1}.session(s).ratemap_even = ratemap_even.session(s).ratemap(counterPlaceCells(s),:);
                        ratemaps.cellN{c,1}.session(s).peakR_even = ratemap_even.session(s).peakR(counterPlaceCells(s),:);
                    end 
                end 
                ratemap_odd.session(s).binPos= vector_binPos_odd;
                ratemap_even.session(s).binPos= vector_binPos_even;
            end 
        end
    end 
    
        if in.whichSessionsToUse(1)==1 && in.whichSessionsToUse(3)==1
            for cell = 1:counterPlaceCells(1)
               [max_val, max_ind] = max(ratemap_odd.session(1).ratemap(cell,:));
               %[max_val, max_ind] = max(ratemap_even.session(1).ratemap(cell,:));
               matrix_peak(cell,1)=cell;
               matrix_peak(cell,2)= max_ind;
            end 

            matrix_order = sortrows(matrix_peak,2);
            
            for cell_i = 1:counterPlaceCells(1)
                cell_num = matrix_order(cell_i,1);
                snakeplot_odd_ses1(cell_i,:) = ratemap_odd.session(1).norm_ratemap(cell_num,:);
                snakeplot_odd_ses2(cell_i,:) = ratemap_odd.session(3).norm_ratemap(cell_num,:);
                snakeplot_even_ses1(cell_i,:) = ratemap_even.session(1).norm_ratemap(cell_num,:);
                snakeplot_even_ses2(cell_i,:) = ratemap_even.session(3).norm_ratemap(cell_num,:);
            end 
            
            ratemaps.session(1).snakeplot_odd = snakeplot_odd_ses1;
            ratemaps.session(2).snakeplot_odd = snakeplot_odd_ses2;
            ratemaps.session(1).snakeplot_even = snakeplot_even_ses1;
            ratemaps.session(2).snakeplot_even = snakeplot_even_ses2;

            odd_trials_ses1 = ceil(length(ratemaps.cellN{1,1}.session(1).Track.trials)/2);
            odd_trials_ses2 = ceil(length(ratemaps.cellN{1,1}.session(3).Track.trials)/2);
            even_trials_ses1 = floor(length(ratemaps.cellN{1,1}.session(1).Track.trials)/2);
            even_trials_ses2 = floor(length(ratemaps.cellN{1,1}.session(3).Track.trials)/2);
            
            if plot_figure == 1
                figure()
                subplot(8,2,[1 3 5])
                s=1;

                imagesc(snakeplot_odd_ses1) %ratemap_odd.session(i).norm_ratemap)
                str_session = sprintf('2IR track ODD - mouse m%d day %d session %d number of trials %d', in.mouse, in.recordingDay,s,odd_trials_ses1);
                ht = title({str_session}, 'FontSize', 8)

                subplot(8,2,7)
                plot(ratemap_odd.session(s).binPos)
                ylabel('binned dwell time (sec)')
                xlim([0 220])

                subplot(8,2,[2 4 6])
                s=3;

                imagesc(snakeplot_odd_ses2) %ratemap_odd.session(i).norm_ratemap)
                str_session = sprintf('2IR track ODD - mouse m%d day %d session %d number of trials %d', in.mouse, in.recordingDay,s,odd_trials_ses2);
                ht = title({str_session}, 'FontSize', 8)

                subplot(8,2,8)
                plot(ratemap_odd.session(s).binPos)
                ylabel('binned dwell time (sec)')
                xlim([0 220])

                subplot(8,2,[9 11 13])
                s=1;

                imagesc(snakeplot_even_ses1) %ratemap_even.session(i).norm_ratemap)
                str_session = sprintf('2IR track EVEN - mouse m%d day %d session %d number of trials %d', in.mouse, in.recordingDay,s,even_trials_ses1);
                ht = title({str_session}, 'FontSize', 8)

                subplot(8,2,15)
                plot(ratemap_even.session(s).binPos)
                ylabel('binned dwell time (sec)')
                xlim([0 220])

                subplot(8,2,[10 12 14])
                s=3;

                imagesc(snakeplot_even_ses2) %ratemap_even.session(i).norm_ratemap)
                str_session = sprintf('2IR track EVEN - mouse m%d day %d session %d number of trials %d', in.mouse, in.recordingDay,s,even_trials_ses2);
                ht = title({str_session}, 'FontSize', 8)

                subplot(8,2,16)
                plot(ratemap_even.session(s).binPos)
                ylabel('binned dwell time (sec)')
                xlim([0 220])
            end 
        
        else 
            missingSession = find (in.whichSessionsToUse(:) == 0);
            if missingSession == 3 % just plot session 1
                ses = 1;
            elseif missingSession == 1 % just plot session 3
                ses = 3;
            end 
             for cell = 1:counterPlaceCells(ses)
               [max_val, max_ind] = max(ratemap_odd.session(ses).ratemap(cell,:));
               matrix_peak(cell,1)=cell;
               matrix_peak(cell,2)= max_ind;
            end 

            matrix_order = sortrows(matrix_peak,2);

            for cell_i = 1:counterPlaceCells(ses)
                cell_num = matrix_order(cell_i,1);
                snakeplot_odd_ses1(cell_i,:) = ratemap_odd.session(ses).norm_ratemap(cell_num,:);
                snakeplot_even_ses1(cell_i,:) = ratemap_even.session(ses).norm_ratemap(cell_num,:);
            end 
            
            if missingSession == 3
                ratemaps.session(1).snakeplot_odd = [];
                ratemaps.session(1).snakeplot_even = [];
                ratemaps.session(2).snakeplot_odd = snakeplot_odd_ses1;
                ratemaps.session(2).snakeplot_even = snakeplot_even_ses1;
            elseif missingSession == 1
                ratemaps.session(2).snakeplot_odd = [];
                ratemaps.session(2).snakeplot_even = [];
                ratemaps.session(1).snakeplot_odd = snakeplot_odd_ses1;
                ratemaps.session(1).snakeplot_even = snakeplot_even_ses1;
            end 
            odd_trials_ses1 = ceil(length(ratemaps.cellN{1,1}.session(ses).Track.trials)/2);
            even_trials_ses1 = floor(length(ratemaps.cellN{1,1}.session(ses).Track.trials)/2);
            
            if plot_figure == 1
                figure()
                subplot(4,2,[1 3 5])
                imagesc(snakeplot_odd_ses1) %ratemap_odd.session(i).norm_ratemap)
                str_session = sprintf('2IR track ODD - mouse m%d day %d session %d number of trials %d', in.mouse, in.recordingDay,ses,odd_trials_ses1);
                ht = title({str_session}, 'FontSize', 8)

                subplot(4,2,7)
                plot(ratemap_odd.session(ses).binPos)
                ylabel('binned dwell time (sec)')
                xlim([0 220])

                subplot(4,2,[2 4 6])
                imagesc(snakeplot_even_ses1) %ratemap_odd.session(i).norm_ratemap)
                str_session = sprintf('2IR track ODD - mouse m%d day %d session %d number of trials %d', in.mouse, in.recordingDay,ses,even_trials_ses1);
                ht = title({str_session}, 'FontSize', 8)

                subplot(4,2,8)
                plot(ratemap_even.session(ses).binPos)
                ylabel('binned dwell time (sec)')
                xlim([0 220])
            end 

        end 
    
    

end 