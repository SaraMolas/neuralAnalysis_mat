% Sara Molas Medina
% 11th October 2023 - updated on January 15th: removed shuffles and check
% all rooms

% function to look at positional shift of place fields in rooms B1 and B2
% in the 2-identical-rooms experiment 

function [posShift] = doPositionalShift (ratemaps,s,in,posShift)
in.room(1).bins = [1:50];
in.room(2).bins = [51:53];
in.room(3).bins = [61:110];
in.room(4).bins = [111:113];
in.room(5).bins = [121:170];
in.room(6).bins = [171:173];
in.room(7).bins = [181:230];
cellsToUse = []; 

for c = 1:length(ratemaps.cellN)
    if ~isfield(ratemaps.cellN{c}.session(s).Track, 'cellOn') || ratemaps.cellN{c}.session(s).Track.cellOn == 0
       continue, 
    end 
    cellsToUse = [cellsToUse, c]; 
end 

rng('shuffle')
randomOrder = randperm(length(cellsToUse)); 
randomCells = cellsToUse(randomOrder); 
for id = 1:length(cellsToUse)
    
    c = cellsToUse(id); 
    randomCell = randomCells(id); 
    posShift.session(s).cell{id}.cellNum = c;
    
    % get place fields in all rooms
    room(1).posField = ratemaps.cellN{c}.session(s).Track.fieldPositions(in.room(1).bins);
    room(2).posField = ratemaps.cellN{c}.session(s).Track.fieldPositions(in.room(3).bins);
    room(3).posField = ratemaps.cellN{c}.session(s).Track.fieldPositions(in.room(5).bins);
    room(4).posField = ratemaps.cellN{c}.session(s).Track.fieldPositions(in.room(7).bins);
    
    % as a control get positional shift in odd vs even trials for all rooms
    ratemapsEven = NaN(1,230); 
    ratemapsOdd = NaN(1,230); 
    counterE = 1;
    counterO = 1;
    for t = 1:length(ratemaps.cellN{c}.session(s).Track.trials)
        if mod(t,2) == 0 %is even
            ratemapsEven(counterE,:) = ratemaps.cellN{c}.session(s).Track.trials{t}.ratemap; 
            counterE = counterE + 1;
        else
            ratemapsOdd(counterO,:) = ratemaps.cellN{c}.session(s).Track.trials{t}.ratemap; 
            counterO = counterO + 1;
        end 
    end 
    % get ratemap of each room for all trials, odd and even
    ratemapOdd = mean(ratemapsOdd, 'omitnan');
    ratemapEven = mean(ratemapsEven, 'omitnan'); 
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
    % get positional shift within room
    for r = 1:4
        if sum(room(r).posField) == 0
            posShift.session(s).cell{id}.room(r) = NaN; 
        else
            [~, peakOdd] = max(room(r).ratemaps.Odd);
            [~, peakEven] = max(room(r).ratemaps.Even);
            posShift.session(s).cell{id}.room(r) = abs(peakOdd - peakEven) * in.spatialBinSize;
        end 
    end 
    % get positional shift across identical rooms
    if (sum(room(2).posField) > 0) && (sum(room(3).posField) > 0)
        [~, peakB1] = max(room(2).ratemaps.All);
        [~, peakB2] = max(room(3).ratemaps.All);
        posShift.session(s).cell{id}.roomsB1B2 = abs(peakB1 - peakB2) * in.spatialBinSize;
    else 
        posShift.session(s).cell{id}.roomsB1B2 = NaN;
    end 
end 

end 