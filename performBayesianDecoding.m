% Sara Molas Medina
% 27 jUNE 2023 - 15th January 2024 (modified so that it uses new room
% coordinates and is more flexible to use)

% Function to decode position on the 2-identical rooms track using bayesian
% decoding - based oon Vandermeer lab Bayesian decoding tutorial

function [BayesianDecoding] = performBayesianDecoding(ratemaps, in, trialData,vr_i,BayesianDecoding,onlyPlaceCells)

% set up some parameters about the rooms
in.rooms{1}.cm.start = 0;
in.rooms{1}.cm.end = 150;
in.rooms{2}.cm.start = 180;
in.rooms{2}.cm.end = 330;
in.rooms{3}.cm.start = 360;
in.rooms{3}.cm.end = 510;
in.rooms{4}.cm.start = 540;
in.rooms{4}.cm.end = 690;
in.rooms{1}.bin.start = 1;
in.rooms{1}.bin.end = 50;
in.rooms{2}.bin.start = 61;
in.rooms{2}.bin.end = 110;
in.rooms{3}.bin.start = 121;
in.rooms{3}.bin.end = 170;
in.rooms{4}.bin.start = 181;
in.rooms{4}.bin.end = 230;
% after calculating decoding errors, we have removed the corridors, so
% rooms have new coordinates
in.rooms{1}.reducBin.start = 1; 
in.rooms{1}.reducBin.end = 50; 
in.rooms{2}.reducBin.start = 51; 
in.rooms{2}.reducBin.end = 100; 
in.rooms{3}.reducBin.start = 101; 
in.rooms{3}.reducBin.end = 150; 
in.rooms{4}.reducBin.start = 151; 
in.rooms{4}.reducBin.end = 200; 

% we are going to do 20% cross-validation - so the data, in this case number of trials, is gonna be split in a 80% portion = training set and a 20% portion = test set, 
% since the first 3 trials tend to be a bit more unstable (according to the Caswell and Romain), I'm going to keep them always in the training set,
% so the test set is going to be a 20% of the trials, from trial 4 to last trial.
if length(ratemaps.cellN{1}.session(vr_i).Track.trials) >= 10
   fractionTest = 0.2;
   TestSize = round(length(ratemaps.cellN{1}.session(vr_i).Track.trials)  * fractionTest);
else
   %fractionTest = 0.2; %round(1/trialData.vr(s).pos.trials(end),1);
   TestSize = 2;
end
TrainSize = length(ratemaps.cellN{1}.session(vr_i).Track.trials)  - TestSize;
%crossval = cvpartition(trialData.vr(s).pos.trials(end), 'Holdout', fractionTest);
BayesianDecodingPosition.session(vr_i).TestSize = TestSize;
BayesianDecodingPosition.session(vr_i).TrainSize = TrainSize;
BayesianDecodingPosition.session(vr_i).seed = rng('shuffle');
BayesianDecodingPosition.session(vr_i).crossValSeq = 3 + randperm(length(ratemaps.cellN{1}.session(vr_i).Track.trials) -3);
    
% create a for-loop to calculate posterior probability for each training set
counterCV = 1;
counterDecoding = 1;
for i = 1:ceil(length(BayesianDecodingPosition.session(vr_i).crossValSeq)/TestSize)
    
    if i == ceil(length(BayesianDecodingPosition.session(vr_i).crossValSeq)/TestSize)
       if i == (length(BayesianDecodingPosition.session(vr_i).crossValSeq)/TestSize)
          testTrials = BayesianDecodingPosition.session(vr_i).crossValSeq(counterCV:(counterCV+TestSize -1));
       else
          remaining = mod(length(BayesianDecodingPosition.session(vr_i).crossValSeq), TestSize);
          testTrials = BayesianDecodingPosition.session(vr_i).crossValSeq(counterCV:(counterCV+remaining -1));
       end          
    else
       testTrials = BayesianDecodingPosition.session(vr_i).crossValSeq(counterCV:(counterCV+TestSize -1));
    end 
        
    %testTrials = BayesianDecodingPosition.session(vr_i).crossValSeq(counterCV:(counterCV+TestSize -1));
    trainTrials = [1:length(ratemaps.cellN{1}.session(vr_i).Track.trials)];
    trainTrials(testTrials) = [];
    counterCV = counterCV + TestSize;
    
    % 1. PREPARE FIRING RATES FOR DECODING 
    %get each cells tuning curve, which for place cells it's a ratemap, but
    % just using 80% of the trials 

    % get place cell ID
    counter=1;
    for c=1:length(ratemaps.cellN)
        if onlyPlaceCells == 1
             if ~isfield(ratemaps.cellN{c}.session(vr_i).Track, 'cellOn') || ratemaps.cellN{c,1}.session(vr_i).Track.cellOn == 0
                continue, 
             end
        end

        % for now I just want to use the place cells in my decoder
        goodCells(counter) = ratemaps.cellN{c,1}.session(vr_i).cellInd(2);
        
        % create ratemap just with the trials of the training set
        counterTrials = 1;
        for t = 1:length(ratemaps.cellN{c,1}.session(vr_i).Track.trials)
            if sum(ismember(trainTrials, t)) == 1
                
%                 if t == trialData.vr(vr_i).pos.trials(end)
%                     completePart = find(~isnan(ratemaps.cellN{c,1}.session(vr_i).Track.trials{t}.ratemap)); %in some trials, session finished on trakc, and then track had nans
%                 else
                    trainingRatemaps(:,1,counterTrials) = ratemaps.cellN{c,1}.session(vr_i).Track.trials{t}.ratemap';
                    counterTrials = counterTrials + 1;
                %end 
            end 
        end 
        % the ratemaps are going to be used as the tuning curves (= encoding model) of these cells
        tuningCurve(:,counter) = mean(trainingRatemaps, 3, 'omitnan');
        counter = counter + 1;
    end 
    BayesianDecodingPosition.session(vr_i).tuningCurve(:,:,counterDecoding) = tuningCurve;
    BayesianDecodingPosition.session(vr_i).set(counterDecoding).Training = trainTrials;
    BayesianDecodingPosition.session(vr_i).set(counterDecoding).Test = testTrials;

    
    % get timeBin index of the test trials
    
    counterDecoding = counterDecoding + 1;
end 


% 1. PREPARE FIRING RATES FOR DECODING 
% Create a matrix where rows is the number of cells and columns is the
% number of time bins. For now just use Place cells
in.timeBin = 1; %0.25; %seconds

spikeTimes = readNPY([in.dataFold filesep 'spike_times.npy']);

% get basic info about spike times
[posData, pos_samp, spikes_session1, spikes_session2, spikes_session3, lastVirmen_session1, lastSleep_session] = concatPos_kilo_SMM(in.fileRoot, spikeTimes,...
    in.times, in.useAllSessions, in.whichSessionsToUse);

% just keep place cell spikes from session 1 and 3 

spikeClusters = readNPY([in.dataFold filesep 'spike_clusters.npy']);
spikesToKeepMask = ismember(spikeClusters,goodCells); % keep only Place cells
PosSampleInRoomVr = find(trialData.vr(vr_i).pos.xy(:,2) == 1 & ((trialData.vr(vr_i).pos.xy(:,1)<in.rooms{1}.cm.end) | (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{2}.cm.start &...
    trialData.vr(vr_i).pos.xy(:,1)<in.rooms{2}.cm.end) | (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{3}.cm.start & trialData.vr(vr_i).pos.xy(:,1)<in.rooms{3}.cm.end) |...
    (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{4}.cm.start & trialData.vr(vr_i).pos.xy(:,1)<in.rooms{4}.cm.end)));
if vr_i == 3
    PosSampleInRoomVr = PosSampleInRoomVr + (lastVirmen_session1 + lastSleep_session);
end 
roomsMask = ismember(pos_samp, PosSampleInRoomVr); % keep only spikes that happened in rooms
if vr_i == 1
    spikes_session = spikes_session1;
elseif vr_i == 3
    spikes_session = spikes_session3;
end 
goodspikes.spikeTimes = double(spikeTimes(roomsMask & spikesToKeepMask & spikes_session))./30000;
goodspikes.pos_sample = pos_samp(roomsMask & spikesToKeepMask & spikes_session);
goodspikes.clusters = spikeClusters(roomsMask & spikesToKeepMask & spikes_session);

% get only positions in the rooms of the track, not in the black
% box or the corridors between rooms. 
trialData.vr(vr_i).pos.track.x = trialData.vr(vr_i).pos.xy(find(trialData.vr(vr_i).pos.xy(:,2) == 1 & ((trialData.vr(vr_i).pos.xy(:,1)<in.rooms{1}.cm.end) |...
    (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{2}.cm.start & trialData.vr(vr_i).pos.xy(:,1)<in.rooms{2}.cm.end) | (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{3}.cm.start &...
    trialData.vr(vr_i).pos.xy(:,1)<in.rooms{3}.cm.end) | (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{4}.cm.start & trialData.vr(vr_i).pos.xy(:,1)<in.rooms{4}.cm.end))),1);
trialData.vr(vr_i).pos.track.ts = trialData.vr(vr_i).pos.ts(find(trialData.vr(vr_i).pos.xy(:,2) == 1 & ((trialData.vr(vr_i).pos.xy(:,1)<in.rooms{1}.cm.end) |...
    (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{2}.cm.start & trialData.vr(vr_i).pos.xy(:,1)<in.rooms{2}.cm.end) | (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{3}.cm.start &...
    trialData.vr(vr_i).pos.xy(:,1)<in.rooms{3}.cm.end) | (trialData.vr(vr_i).pos.xy(:,1)>in.rooms{4}.cm.start & trialData.vr(vr_i).pos.xy(:,1)<in.rooms{4}.cm.end))),1);

% when sampling in time windows how do I deal with the fact that I just
% want to look at the rooms, not at the black box or the corridors? 
binEdges = trialData.vr(vr_i).pos.track.ts(1):in.timeBin:trialData.vr(vr_i).pos.track.ts(end);
binCenters = binEdges(1:end-1)+in.timeBin/2;
trialsBinned = interp1(trialData.vr(vr_i).pos.ts, trialData.vr(vr_i).pos.trials,binCenters);
% create a matrix where rows = number of cells and columns = number of time
% bins. in each bin you get spike count. 
for c = 1:length(goodCells)
    placeCell = goodCells(c); 
    
    % get the spikes for that cell 
    spikeTS = goodspikes.spikeTimes(goodspikes.clusters == placeCell);
    cellTimeMatrix(c, :)=histc(spikeTS, binEdges); 
    cellTimeMatrix(c, end-1) = cellTimeMatrix(c, end-1)+cellTimeMatrix(c, end); %vandermeer does this bc of last bin of histc() gotcha?
end 
cellTimeMatrix = cellTimeMatrix(:,1: end-1);

% now we have a numberCells*numberTimeBins matrix - let's look at it
% figure()
% imagesc(binCenters,1:length(goodCells),cellTimeMatrix)
% set(gca,'FontSize',16); xlabel('time(s)'); ylabel('cell #');

% 2 - CREATE OUR PRIOR PROBABILITY (P(X))
% for now our prior will be a uniform occupation matrix, in the future try
% it with the real occupancy dwell times
nBins = size(tuningCurve,1);
nCells = size(tuningCurve,2);
occUniform = repmat (1/nBins, [nBins 1]);
%dwellTime = ratemaps.vr(vr_i).Track.posBinData([1:50, 57:107,114:164,171:220])';
%pDwellTime = dwellTime ./(sum(dwellTime));

% 3 - RUN THE BAYESIAN DECODER - p(x|n) = p(n|x) * p(x) / p(n)
% iterate through the training sets and store each posterior 
counterDecoding = 1;
BayesianDecodingPosition.session(vr_i).posteriorTest = NaN(nBins, length(binCenters));
for i = 1:ceil(length(BayesianDecodingPosition.session(vr_i).crossValSeq)/TestSize)
    len = length(binCenters);
    p = NaN(length(binCenters),nBins);
    nActiveNeurons = sum(cellTimeMatrix > 0);
    tuningCurve = BayesianDecodingPosition.session(vr_i).tuningCurve(:,:,counterDecoding);
    for iB = 1:nBins
        tempProd = nansum(log(repmat(tuningCurve(iB,:)',1,len) .^ cellTimeMatrix));
        tempSum = exp(-in.timeBin * nansum(tuningCurve(iB,:),2));
        p(:,iB) = exp(tempProd)*tempSum*occUniform(iB);
    end 
    BayesianDecodingPosition.session(vr_i).posterior(:,:,counterDecoding) = p;
    BayesianDecodingPosition.session(vr_i).posteriorNormalized(:,:,counterDecoding) = p./repmat(sum(p,2),1,nBins); % renormalized to 1 total probability 
    p(nActiveNeurons < 1,:) = 0;
    BayesianDecodingPosition.session(vr_i).posteriorIgnoringInactiveBins(:,:,counterDecoding) = p;
    
    testIdx = [];
    for t = 1:length(BayesianDecodingPosition.session(vr_i).set(counterDecoding).Test)
        testTrials = BayesianDecodingPosition.session(vr_i).set(counterDecoding).Test;
        trial = testTrials(t); 
        trialIdx = find(trialsBinned == trial);
        testIdx = [testIdx, trialIdx];
    end 
    BayesianDecodingPosition.session(vr_i).posteriorTestSet(counterDecoding).p(:,:) = p(testIdx,:)';
    BayesianDecodingPosition.session(vr_i).testBinsIndex(counterDecoding).testIdx(:,:) = testIdx;
    BayesianDecodingPosition.session(vr_i).posteriorTest(:,testIdx) = p(testIdx,:)';
    counterDecoding = counterDecoding + 1;
end 



% 4 - VISUALIZE THE RESULTS

% Plot where the x-axis is time, the y-axis is position. In imagesc, plot the probability of the mouse's positon, and then with a dot? or line? in plot, plot the
% real position of the mouse
%xBinned = interp1(trialData.vr(1).track.pos.ts, trialData.vr(1).track.pos.xy(:,1),binCenters);
xBinned = interp1(trialData.vr(vr_i).pos.ts, trialData.vr(vr_i).pos.xy(:,1),binCenters, 'nearest');
yBinned = interp1(trialData.vr(vr_i).pos.ts, trialData.vr(vr_i).pos.xy(:,2),binCenters, 'nearest');
speedBinned = interp1(trialData.vr(vr_i).pos.ts, trialData.vr(vr_i).pos.speed,binCenters, 'nearest');
notRoomsPos = find(yBinned == 2 | ((xBinned > in.rooms{1}.cm.end & xBinned <in.rooms{2}.cm.start)|(xBinned > in.rooms{2}.cm.end & xBinned <in.rooms{3}.cm.start)|...
    (xBinned > in.rooms{3}.cm.end & xBinned <in.rooms{4}.cm.start)|(xBinned>in.rooms{4}.cm.end)));
xBinned(notRoomsPos) = NaN;
stationaryPos = find(speedBinned < 2);
xBinned(stationaryPos) = NaN;
earlyTrials = find(trialsBinned <= 3);
xBinned(earlyTrials) = NaN;

xBinnedBins = round(xBinned ./ 3); %discretize(xBinned, size(p,1));
% need to modify the x-binned positions so I only take into account the
% rooms when calculating decoding error
roomB1Pos = find(xBinnedBins>=in.rooms{2}.bin.start & xBinnedBins<= in.rooms{2}.bin.end);%[1:50, 57:107,114:164,171:220]
xBinnedBins(roomB1Pos) = xBinnedBins(roomB1Pos) - 10;%6;
roomB2Pos = find(xBinnedBins>=in.rooms{3}.bin.start & xBinnedBins<= in.rooms{3}.bin.end);
xBinnedBins(roomB2Pos) = xBinnedBins(roomB2Pos) - 20;%12;
roomCPos = find(xBinnedBins>=in.rooms{4}.bin.start & xBinnedBins<= in.rooms{4}.bin.end);
xBinnedBins(roomCPos) = xBinnedBins(roomCPos) - 30; %18;
% get decoded position 
decodedPosition = NaN(size(BayesianDecodingPosition.session(vr_i).posteriorTest,2),3);
for i = 1:size(BayesianDecodingPosition.session(vr_i).posteriorTest,2)
    % find decoded position based on highest p value
    [maxP, pos] = max(BayesianDecodingPosition.session(vr_i).posteriorTest(:,i)); % this should be  max(BayesianDecodingPosition.session(vr_i).posteriorTest(:,i)) no? 
    if isnan(maxP)
        pos = NaN;
%         if any(ismember(notRoomsPos, i)) | any(ismember(stationaryPos, i)) | any(ismember(earlyTrials, i))
%            
%         else
%             disp('p is NaN')
%         end 
    end 
    decodedPosition(i,1) = pos; 
    decodedPosition(i,2) = binCenters(i);
    % find decoded position based on highest p value inside the room with
    % highest total probability
    prob = NaN(4,1);
    for r = 1:4
        prob(r) = sum(BayesianDecodingPosition.session(vr_i).posteriorTest([in.rooms{r}.bin.start:in.rooms{r}.bin.end],i),'omitnan');
    end 
    [~, room] = max(prob); 
    [~, pos] = max(BayesianDecodingPosition.session(vr_i).posteriorTest([in.rooms{room}.bin.start:in.rooms{room}.bin.end],i));
    decodedPosition(i,3) = in.rooms{room}.bin.start + pos - 1; 
end 

decodedPosition(notRoomsPos,[1,3]) = NaN;
decodedPosition(stationaryPos,[1,3]) = NaN;
decodedPosition(earlyTrials,[1,3]) = NaN;
% decodedPosition(notRoomsPos,3) = NaN;
% decodedPosition(stationaryPos,3) = NaN;
% decodedPosition(earlyTrials,3) = NaN;
roomB1Pos = find(decodedPosition(:,1)>=in.rooms{2}.bin.start & decodedPosition(:,1)<= in.rooms{2}.bin.end);%[1:50, 57:107,114:164,171:220]
decodedPosition(roomB1Pos,1) = decodedPosition(roomB1Pos,1) - 10;%6;
roomB2Pos = find(decodedPosition(:,1)>=in.rooms{3}.bin.start & decodedPosition(:,1)<= in.rooms{3}.bin.end);
decodedPosition(roomB2Pos,1) = decodedPosition(roomB2Pos,1) - 20;%12;
roomCPos = find(decodedPosition(:,1)>=in.rooms{4}.bin.start & decodedPosition(:,1)<= in.rooms{4}.bin.end);
decodedPosition(roomCPos,1) = decodedPosition(roomCPos,1) - 30; %18;

roomB1Pos = find(decodedPosition(:,3)>=in.rooms{2}.bin.start & decodedPosition(:,3)<= in.rooms{2}.bin.end);%[1:50, 57:107,114:164,171:220]
decodedPosition(roomB1Pos,3) = decodedPosition(roomB1Pos,3) - 10;%6;
roomB2Pos = find(decodedPosition(:,3)>=in.rooms{3}.bin.start & decodedPosition(:,3)<= in.rooms{3}.bin.end);
decodedPosition(roomB2Pos,3) = decodedPosition(roomB2Pos,3) - 20;%12;
roomCPos = find(decodedPosition(:,3)>=in.rooms{4}.bin.start & decodedPosition(:,3)<= in.rooms{4}.bin.end);
decodedPosition(roomCPos,3) = decodedPosition(roomCPos,3) - 30; %18;

% get decoding error as if track was circular
%Plot decoding error across rooms
BayesianDecoding.session(vr_i).wholeTrack.decErrorMaxP.Cm = NaN(length(xBinned),1);
BayesianDecoding.session(vr_i).wholeTrack.decErrorMaxP.Sqrt = NaN(length(xBinned),1);
BayesianDecoding.session(vr_i).wholeTrack.decErrorSumP.Cm = NaN(length(xBinned),1);
BayesianDecoding.session(vr_i).wholeTrack.decErrorSumP.Sqrt = NaN(length(xBinned),1);
for r = 1:4
    BayesianDecoding.session(vr_i).rooms{r}.decErrorMaxP.Sqrt = NaN(length(xBinned),1);
    BayesianDecoding.session(vr_i).rooms{r}.decErrorMaxP.Cm = NaN(length(xBinned),1);
    BayesianDecoding.session(vr_i).rooms{r}.decErrorSumP.Sqrt = NaN(length(xBinned),1);
    BayesianDecoding.session(vr_i).rooms{r}.decErrorSumP.Cm = NaN(length(xBinned),1);
end 
for x = 1:length(xBinnedBins)
    % based on position decoded with max P value
    if abs(xBinnedBins(x)-decodedPosition(x,1)) > (in.rooms{4}.bin.end/ 2) % check that xBinned is in cm and decoded position too
        BayesianDecoding.session(vr_i).wholeTrack.decErrorMaxP.Cm(x) =abs(abs(xBinnedBins(x)-decodedPosition(x,1))-in.rooms{4}.bin.end);
        BayesianDecoding.session(vr_i).wholeTrack.decErrorMaxP.Sqrt(x) = sqrt((abs(abs(xBinnedBins(x)-decodedPosition(x,1))-in.rooms{4}.bin.end)).^2);
    else
        BayesianDecoding.session(vr_i).wholeTrack.decErrorMaxP.Cm(x) = abs(xBinnedBins(x)-decodedPosition(x,1));
        BayesianDecoding.session(vr_i).wholeTrack.decErrorMaxP.Sqrt(x) = sqrt(abs(xBinnedBins(x)-decodedPosition(x,1)) .^2);
    end 
    for r = 1:4
        if xBinned(x) >= in.rooms{r}.cm.start && xBinned(x) <= in.rooms{r}.cm.end % if mouse in room A
            BayesianDecoding.session(vr_i).rooms{r}.decErrorMaxP.Sqrt(x) = BayesianDecoding.session(vr_i).wholeTrack.decErrorMaxP.Sqrt(x);
            BayesianDecoding.session(vr_i).rooms{r}.decErrorMaxP.Cm(x) = BayesianDecoding.session(vr_i).wholeTrack.decErrorMaxP.Cm(x);
        end 
    end 
    % based on position decoded with max P value
    if abs(xBinnedBins(x)-decodedPosition(x,3)) > (in.rooms{4}.bin.end/ 2) % check that xBinned is in cm and decoded position too
        BayesianDecoding.session(vr_i).wholeTrack.decErrorSumP.Cm(x) =abs(abs(xBinnedBins(x)-decodedPosition(x,3))-in.rooms{4}.bin.end);
        BayesianDecoding.session(vr_i).wholeTrack.decErrorSumP.Sqrt(x) = sqrt((abs(abs(xBinnedBins(x)-decodedPosition(x,3))-in.rooms{4}.bin.end)).^2);
    else
        BayesianDecoding.session(vr_i).wholeTrack.decErroSumP.Cm(x) = abs(xBinnedBins(x)-decodedPosition(x,3));
        BayesianDecoding.session(vr_i).wholeTrack.decErrorSumP.Sqrt(x) = sqrt(abs(xBinnedBins(x)-decodedPosition(x,3)) .^2);
    end 
    for r = 1:4
        if xBinned(x) >= in.rooms{r}.cm.start && xBinned(x) <= in.rooms{r}.cm.end % if mouse in room A
            BayesianDecoding.session(vr_i).rooms{r}.decErrorSumP.Sqrt(x) = BayesianDecoding.session(vr_i).wholeTrack.decErrorSumP.Sqrt(x);
            BayesianDecoding.session(vr_i).rooms{r}.decErrorSumP.Cm(x) = BayesianDecoding.session(vr_i).wholeTrack.decErrorSumP.Cm(x);
        end 
    end 
end 

figure()
im = imagesc(p', [0 1]);
xlabel ('Time (s)');
ylabel ('Position along the track');
yticks ([25, 75, 125, 175]);
yticklabels ({'Room A', 'Room B1', 'Room B2', 'Room C'});
title_str = sprintf('Actual(green) vs decoded(red) position in mouse m%d day %d session %d', in.mouse, in.recordingDay, vr_i);
title(title_str);
hold on 
plot(1:length(decodedPosition(:,2)), decodedPosition(:,1), 'o', 'MarkerFaceColor', 'red');
plot(1:length(decodedPosition(:,2)), xBinnedBins(:), 'o',  'MarkerFaceColor', 'green');
hold off
set(gca, 'FontSize', 30)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
figname = sprintf('Analysis\BayesDecPos_%d.png', vr_i);
saveas(gcf,[in.filePath filesep figname]);
close all


% Check what is the actual room and then what room is being decoded based
% on highest probability value in the whole track and highest room probability (calculated by adding up all the probabilities in each room)
roomVectors.actual.room = NaN(length(xBinned),1);
roomVectors.decoded.maxP.room = NaN(length(xBinned),1);
roomVectors.decoded.sumP.room = NaN(length(xBinned),1);
roomVectors.actual.location = NaN(length(xBinned),1);
roomVectors.decoded.maxP.location = NaN(length(xBinned),1);
roomVectors.decoded.sumP.location = NaN(length(xBinned),1);

for x = 1:length(xBinned)
    for r = 1:4 
        if xBinned(x) >= in.rooms{r}.cm.start && xBinned(x) <= in.rooms{r}.cm.end  % if mouse in room A
            roomVectors.actual.room(x)= r;
        end 
    end  
    if (xBinned(x) >= in.rooms{1}.cm.start  && xBinned(x) <= (in.rooms{1}.cm.start+50)) || (xBinned(x) >= in.rooms{2}.cm.start  && xBinned(x) <= (in.rooms{2}.cm.start+50)) ||...
            (xBinned(x) >= in.rooms{3}.cm.start  && xBinned(x) <= (in.rooms{3}.cm.start+50)) || (xBinned(x) >= in.rooms{4}.cm.start  && xBinned(x) <= (in.rooms{4}.cm.start+50))
        roomVectors.actual.location(x)=1; % start
    elseif (xBinned(x) >= (in.rooms{1}.cm.start+50) && xBinned(x) <= (in.rooms{1}.cm.start+100)) || (xBinned(x) >= (in.rooms{2}.cm.start+50) && xBinned(x) <= (in.rooms{2}.cm.start+100)) ||...
            (xBinned(x) >= (in.rooms{3}.cm.start+50) && xBinned(x) <= (in.rooms{3}.cm.start+100))  || (xBinned(x) >= (in.rooms{4}.cm.start+50) && xBinned(x) <= (in.rooms{4}.cm.start+100)) 
        roomVectors.actual.location(x)=2; % middle
    elseif (xBinned(x) >= (in.rooms{1}.cm.start+100) && xBinned(x) <= (in.rooms{1}.cm.end))  || (xBinned(x) >= (in.rooms{2}.cm.start+100) && xBinned(x) <= (in.rooms{2}.cm.end)) ||...
            (xBinned(x) >= (in.rooms{3}.cm.start+100) && xBinned(x) <= (in.rooms{3}.cm.end)) || (xBinned(x) >= (in.rooms{4}.cm.start+100) && xBinned(x) <= (in.rooms{4}.cm.end))
        roomVectors.actual.location(x)=3; % end
    end 
    for r = 1:4
        if decodedPosition(x,1)>= in.rooms{r}.reducBin.start  && decodedPosition(x,1)<=in.rooms{r}.reducBin.end % if max P value says mouse in room R
            roomVectors.decoded.maxP.room(x)=r;
        end 
    end 
    if (decodedPosition(x,1)>=in.rooms{1}.reducBin.start && decodedPosition(x,1)<= (in.rooms{1}.reducBin.start + 50 /3)) ||...
            (decodedPosition(x,1)>in.rooms{2}.reducBin.start && decodedPosition(x,1)<= (in.rooms{2}.reducBin.start + 50 /3)) ||...
            (decodedPosition(x,1)>in.rooms{3}.reducBin.start && decodedPosition(x,1)<= (in.rooms{3}.reducBin.start + 50 /3)) ||...
            (decodedPosition(x,1)>in.rooms{4}.reducBin.start && decodedPosition(x,1)<=(in.rooms{4}.reducBin.start + 50 /3))% if max P value says mouse in start of room
       roomVectors.decoded.maxP.location(x)=1;
    elseif (decodedPosition(x,1)>=(in.rooms{1}.reducBin.start+(50 /3)) && decodedPosition(x,1)<= (in.rooms{1}.reducBin.start + 2*(50 /3)) ||...
            (decodedPosition(x,1)>(in.rooms{2}.reducBin.start+(50 /3)) && decodedPosition(x,1)<= (in.rooms{2}.reducBin.start + 2*(50 /3))) ||...
            (decodedPosition(x,1)>(in.rooms{3}.reducBin.start+(50 /3)) && decodedPosition(x,1)<= (in.rooms{3}.reducBin.start + 2*(50 /3))) ||...
            (decodedPosition(x,1)>(in.rooms{4}.reducBin.start+(50 /3)) && decodedPosition(x,1)<=(in.rooms{4}.reducBin.start + 2*(50 /3)))) 
        roomVectors.decoded.maxP.location(x)=2;
    elseif (decodedPosition(x,1)>=(in.rooms{1}.reducBin.start+2*(50 /3)) && decodedPosition(x,1)<= in.rooms{1}.reducBin.end) ||...
            (decodedPosition(x,1)>(in.rooms{2}.reducBin.start+2*(50 /3)) && decodedPosition(x,1)<= in.rooms{2}.reducBin.end) ||...
            (decodedPosition(x,1)>(in.rooms{3}.reducBin.start+2*(50 /3)) && decodedPosition(x,1)<= in.rooms{3}.reducBin.end) ||...
            (decodedPosition(x,1)>(in.rooms{4}.reducBin.start+2*(50 /3)) && decodedPosition(x,1)<=in.rooms{4}.reducBin.end)% if max P value says mouse in end of room
        roomVectors.decoded.maxP.location(x)=3;
    end 
end 

for x=1:length(xBinned)
    room = NaN(4,1); 
    for r = 1:4
        room(r) = sum(p(x,[in.rooms{r}.reducBin.start:in.rooms{r}.reducBin.end]));
    end 
    [maxp, maxRoom] = max(room);
    roomVectors.decoded.sumP.room(x)=maxRoom;
    location = NaN(3,1);
    location(1) = sum(p(x,[in.rooms{maxRoom}.reducBin.start:in.rooms{maxRoom}.reducBin.start+ceil(50/3)]));
    location(2) = sum(p(x,[in.rooms{maxRoom}.reducBin.start+ceil(50/3):in.rooms{maxRoom}.reducBin.start+2*ceil(50/3)]));
    location(3) = sum(p(x,[in.rooms{maxRoom}.reducBin.start+2*ceil(50/3):in.rooms{maxRoom}.reducBin.end]));
    [maxp, maxLocation] = max(location);
    roomVectors.decoded.sumP.location(x)=maxLocation;
end 

if max(xBinnedBins) > in.rooms{4}.reducBin.end
    disp('xBinned bin is larger than 200')
end
BayesianDecodingPosition.session(vr_i).roomVectors = roomVectors;
BayesianDecodingPosition.session(vr_i).xbinnedBins = xBinnedBins;
BayesianDecodingPosition.session(vr_i).decodedPosition = decodedPosition;
BayesianDecoding.session(vr_i).results = BayesianDecodingPosition.session(vr_i); 

end 
