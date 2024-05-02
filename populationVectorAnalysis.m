% Sara Molas Medina
% 27 June 2023
% updated on 1st November 2023 so it will work with ratemaps that have 1
% time bin per corridor

function [populationVectorAnalysis] = populationVectorAnalysis (ratemaps, in, trialData)

    
for s=1:2:3
    if in.whichSessionsToUse(s) == 0
       continue, 
    end 
        
            %perform 10-fold cross validation
    % take 10% of the trials, ignoring the first 3 trials of a session (since it's still not that stable there, so they should remain in the training set) and then do this iteratively
    %numberTrialsCV = trialData.vr(s).pos.trials(end) - 4;
    if length(ratemaps.cellN{1}.session(s).Track.trials) >= 10
        fractionTest = 0.2;
        TestSize = round(length(ratemaps.cellN{1}.session(s).Track.trials) * fractionTest);
    else
         %round(1/trialData.vr(s).pos.trials(end),1);
        TestSize = 2;
    end
    
    TrainSize = length(ratemaps.cellN{1}.session(s).Track.trials) - TestSize;
    %crossval = cvpartition(trialData.vr(s).pos.trials(end), 'Holdout', fractionTest);
    populationVectorAnalysis.session(s).TestSize = TestSize;
    populationVectorAnalysis.session(s).TrainSize = TrainSize;
    populationVectorAnalysis.session(s).seed = rng('shuffle');
    populationVectorAnalysis.session(s).crossValSeq = 3 + randperm(length(ratemaps.cellN{1}.session(s).Track.trials)-3);
    
    % loop over cross-validation training and test sets
    counterCV = 1;
    counterCrossVal = 1;
    maxCounterCrossVal = ceil(length(populationVectorAnalysis.session(s).crossValSeq)/TestSize);
    populationVectorAnalysis.session(s).ratemaps.TrainTestCorrelation([1:203],[1:203],[1:maxCounterCrossVal]) = NaN;
    populationVectorAnalysis.session(s).ratemaps.TrainTestCorrelationHigh([1:203],[1:203],[1:maxCounterCrossVal]) = NaN;
    populationVectorAnalysis.session(s).ratemaps.TrainTestCorrelationLow([1:203],[1:203],[1:maxCounterCrossVal]) = NaN;
    for i = 1:ceil(length(populationVectorAnalysis.session(s).crossValSeq)/TestSize)
        if i == ceil(length(populationVectorAnalysis.session(s).crossValSeq)/TestSize)
            if i == (length(populationVectorAnalysis.session(s).crossValSeq)/TestSize)
                testTrials = populationVectorAnalysis.session(s).crossValSeq(counterCV:(counterCV+TestSize -1));
            else
                remaining = mod(length(populationVectorAnalysis.session(s).crossValSeq), TestSize);
                testTrials = populationVectorAnalysis.session(s).crossValSeq(counterCV:(counterCV+remaining -1));
            end 
                
        else
            testTrials = populationVectorAnalysis.session(s).crossValSeq(counterCV:(counterCV+TestSize -1));
        end 
        counterCV = counterCV + TestSize;
        
        trainTrials = [1:length(ratemaps.cellN{1}.session(s).Track.trials)];
        trainTrials(testTrials) = [];
        populationVectorAnalysis.session(s).sets(counterCrossVal).Training = trainTrials;
        populationVectorAnalysis.session(s).sets(counterCrossVal).Test = testTrials;
                
        counterCell = 1;
        counterCellLow = 1;
        counterCellHigh = 1;
        
        % first get the ratemaps of the place cells active in that session
        for c = 1:length(ratemaps.cellN)
            if s==3 &&  length(ratemaps.cellN{c,1}.session)==1
                continue, 
            end 
            if ~isfield(ratemaps.cellN{c}.session(s).Track, 'cellOn') || ratemaps.cellN{c,1}.session(s).Track.cellOn == 0 % only take place cells
                continue,
            end 
            mRatemap_A = NaN(1, length(ratemaps.cellN{c,1}.session(s).Track.ratemap));
            mRatemap_B = NaN(1, length(ratemaps.cellN{c,1}.session(s).Track.ratemap));
            counter_A = 1;
            counter_B = 1;
            for trial_i = 1: length(ratemaps.cellN{c}.session(s).Track.trials)
                if sum(ismember(trainTrials,trial_i)) == 1
                   mRatemap_A(counter_A,:) = ratemaps.cellN{c}.session(s).Track.trials{trial_i}.ratemap;
                   counter_A = counter_A + 1;
                elseif sum(ismember(testTrials,trial_i)) == 1
                    mRatemap_B(counter_B,:) = ratemaps.cellN{c}.session(s).Track.trials{trial_i}.ratemap;
                    counter_B = counter_B + 1;
                else 
                    disp('Problem with trials in training and test sets - Population Vector Analysis')
                end 
            end 
            populationVectorAnalysis.session(s).ratemaps.TrainTrials(counterCell,:) = (mean(mRatemap_A, 'omitnan')) ./ max(mean(mRatemap_A, 'omitnan'));
            populationVectorAnalysis.session(s).ratemaps.TestTrials(counterCell,:) = (mean(mRatemap_B, 'omitnan')) ./ max(mean(mRatemap_B, 'omitnan'));
            counterCell = counterCell + 1;
        end 

        NaNCol = [];
        for col = 1:size(mRatemap_B, 2)
            if isnan(populationVectorAnalysis.session(s).ratemaps.TrainTrials(:,col))
                NaNCol = [NaNCol, col];
            end 
            if isnan(populationVectorAnalysis.session(s).ratemaps.TestTrials(:,col))
                NaNCol = [NaNCol, col];
            end 
        end 
        colsMask = [1:size(mRatemap_B, 2)];
        if ~isempty(NaNCol)
            colsMask(NaNCol)= [];
        end 
        
        populationVectorAnalysis.session(s).ratemaps.TrainTestCorrelation(colsMask,colsMask,counterCrossVal) = corr (populationVectorAnalysis.session(s).ratemaps.TrainTrials(:,colsMask), ...
          populationVectorAnalysis.session(s).ratemaps.TestTrials(:,colsMask), 'rows', 'complete'); 
      populationVectorAnalysis.session(s).ratemaps.colsMask{counterCrossVal} = colsMask;
       counterCrossVal = counterCrossVal +1;
    end 
    
        populationVectorAnalysis.session(s).ratemaps.meanTrainTestCorrelation = mean(populationVectorAnalysis.session(s).ratemaps.TrainTestCorrelation,3, 'omitnan');
        % January 14th: added 'omitnan' when computing the mean
end        

end 
