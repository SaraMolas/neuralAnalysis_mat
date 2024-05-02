% Sara Molas Medina
% 11th October 2023 - modified 14th January 2024 to include all rooms

% function to look at rate remapping between rooms B1 and B2 in
% 2-identical-rooms experiment 

function [rateRemap] = doRateRemapping (ratemaps,s,rateRemap, in, useCellsWithOnePFperRoom)
% useCellsWithOnePFperRoom == 0, use cells with one pf in room B1 or B2, if
% == 1, use cells with two pfs, one in room B1 and one in room B2
trials_A = [1:2:length(ratemaps.vr(s).Track.trial)];
trials_B = [2:2:length(ratemaps.vr(s).Track.trial)];

for c = 1:length(ratemaps.cellN)

    % for now skip if it's not a PC
    if ~isfield(ratemaps.cellN{c}.session(s).Track, 'cellOn') || ratemaps.cellN{c}.session(s).Track.cellOn == 0
        continue, 
    end 
    
    mRatemap_A = NaN(ceil(length(ratemaps.cellN{c}.session(s).Track.trials)/2), size(ratemaps.cellN{c}.session(s).Track.ratemap,2));
    mRatemap_B = NaN(ceil(length(ratemaps.cellN{c}.session(s).Track.trials)/2), size(ratemaps.cellN{c}.session(s).Track.ratemap,2));
    counter_A = 1;
    counter_B = 1;
    for trial_i = 1: length(ratemaps.cellN{c}.session(s).Track.trials) %trialData.vr(s).pos.trials(end)
        if ismember(trial_i, trials_A)
            mRatemap_A(counter_A,:) = ratemaps.cellN{c}.session(s).Track.trials{trial_i}.ratemap;
            counter_A = counter_A + 1;
        elseif ismember(trial_i, trials_B)
            mRatemap_B(counter_B,:) = ratemaps.cellN{c}.session(s).Track.trials{trial_i}.ratemap;
            counter_B = counter_B + 1;
        end 
    end 
    ratemap_A = mean(mRatemap_A, 'omitnan');
    ratemap_B = mean(mRatemap_B, 'omitnan');  
    
    ratemapOdd = ratemap_A; % ratemaps.cellN{c}.session(s).Track.ratemapOdd
    ratemapEven = ratemap_B;%ratemaps.cellN{c}.session(s).Track.ratemapEven
    
    % if only pf in one room but not the other, take mean firing rate in the exisiting pf as the difference
    room(1).ratemaps.Odd = ratemapOdd(in.room(1).bins);
    room(1).ratemaps.Even = ratemapEven(in.room(1).bins);
    room(1).ratemaps.All = ratemaps.cellN{c}.session(s).Track.ratemap(in.room(1).bins);
    room(2).ratemaps.Odd = ratemapOdd(in.room(3).bins);
    room(2).ratemaps.Even = ratemapEven(in.room(3).bins);
    room(2).ratemaps.All = ratemaps.cellN{c}.session(s).Track.ratemap(in.room(3).bins);
    room(3).ratemaps.Odd = ratemapOdd(in.room(5).bins); 
    room(3).ratemaps.Even = ratemapEven(in.room(5).bins);
    room(3).ratemaps.All = ratemaps.cellN{c}.session(s).Track.ratemap(in.room(5).bins);
    room(4).ratemaps.Odd = ratemapOdd(in.room(7).bins);
    room(4).ratemaps.Even = ratemapEven(in.room(7).bins);
    room(4).ratemaps.All = ratemaps.cellN{c}.session(s).Track.ratemap(in.room(7).bins);
    room(1).posField = ratemaps.cellN{c}.session(s).Track.fieldPositions(in.room(1).bins);
    room(2).posField = ratemaps.cellN{c}.session(s).Track.fieldPositions(in.room(3).bins);
    room(3).posField = ratemaps.cellN{c}.session(s).Track.fieldPositions(in.room(5).bins);
    room(4).posField = ratemaps.cellN{c}.session(s).Track.fieldPositions(in.room(7).bins);
    for r = 1:4
        if sum(room(r).posField) == 0
            room(r).meanOdd = NaN;
            room(r).meanEven = NaN;
            room(r).maxOdd = NaN;
            room(r).maxEven = NaN;
            room(r).mean = NaN; 
            room(r).max = NaN;
        else 
            room(r).meanOdd = mean(room(r).ratemaps.Odd(logical(room(r).posField)), 'omitnan');
            room(r).meanEven = mean(room(r).ratemaps.Even(logical(room(r).posField)), 'omitnan');
            room(r).maxOdd = max(room(r).ratemaps.Odd(logical(room(r).posField)));
            room(r).maxEven = max(room(r).ratemaps.Even(logical(room(r).posField)));
            room(r).mean = mean(room(r).ratemaps.All(logical(room(r).posField)),'omitnan');
            room(r).max = max(room(r).ratemaps.All(logical(room(r).posField)));
        end
    end 
    % get rate remapping within room
    for r = 1:4 
        if sum(room(r).posField) == 0
            rateRemap.session(s).cellN{c}.room(r).mean = NaN;
            rateRemap.session(s).cellN{c}.room(r).max = NaN;
        else
            rateRemap.session(s).cellN{c}.room(r).mean = abs(room(r).meanOdd - room(r).meanEven)/(room(r).meanOdd + room(r).meanEven);
            rateRemap.session(s).cellN{c}.room(r).max = abs(room(r).maxOdd - room(r).maxEven)/(room(r).maxOdd + room(r).maxEven);
        end
    end
    % get rate remapping across identical rooms
    if (sum(room(2).posField) > 0) && (sum(room(3).posField) > 0)
        rateRemap.session(s).cellN{c}.roomsB1B2.mean = abs(room(2).mean - room(3).mean)/(room(2).mean + room(3).mean);
        rateRemap.session(s).cellN{c}.roomsB1B2.max = abs(room(2).max - room(3).max)/(room(2).max + room(3).max);
    else
        rateRemap.session(s).cellN{c}.roomsB1B2.mean = NaN;
        rateRemap.session(s).cellN{c}.roomsB1B2.max = NaN;
        if (sum(room(2).posField) > 0)
            rateRemap.session(s).cellN{c}.roomsB1.mean = abs(room(2).mean - 0)/(room(2).mean + 0);
            rateRemap.session(s).cellN{c}.roomsB1.max = abs(room(2).max - 0)/(room(2).max + 0);
        elseif (sum(room(3).posField) > 0)
            rateRemap.session(s).cellN{c}.roomsB2.mean = abs(room(3).mean - 0)/(room(3).mean + 0);
            rateRemap.session(s).cellN{c}.roomsB2.max = abs(room(3).max - 0)/(room(3).max + 0);
        end 
    end 
end
end 