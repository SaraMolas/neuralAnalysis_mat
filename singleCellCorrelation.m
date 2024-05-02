% Sara Molas Medina
%14th September 2023
% updated on 31st October 2023 - so it works for ratemaps with time-binned
% corridors. also no need to compute even and odd ratemaps bc that's been
% done before

% Function to compute single cell correlations in 2-identical rooms task 
% just wanted to make it nicer than the compute intrasession_ocrrelations
% in the main runSaraAnalysis.m master script

function [scCor, in, ratemaps] = singleCellCorrelation (scCor, trialData,ratemaps,s, in)

% STEP 1: Compute place cell stability for each place cell and overall for
% that session


cellStability = [];
counter = 1;
ratemapsOdd = NaN(1,230);
ratemapsEven = NaN(1,230);
for c = 1:length(ratemaps.cellN)
    if s== 3 && length(ratemaps.cellN{c,1}.session)
        ratemaps.cellN{c}.session(s).Track.cellStability = NaN;
        continue, 
    end 
    if ratemaps.cellN{c,1}.session(s).Track.cellOn == 0
        ratemaps.cellN{c}.session(s).Track.cellStability = NaN;
        continue, 
    end 
    
    [ratemaps.cellN{c}.session(s).Track.cellStability, ~] = spatialCorrelation(ratemaps.cellN{c,1}.session(s).Track.ratemapEven,...
        ratemaps.cellN{c,1}.session(s).Track.ratemapOdd,false);
    ratemapsOdd(counter,:) = ratemaps.cellN{c,1}.session(s).Track.ratemapOdd;
    ratemapsEven(counter,:) = ratemaps.cellN{c,1}.session(s).Track.ratemapEven;
    cellIDs(counter) = c; 
    counter = counter+1; 
    if ratemaps.cellN{c}.session(s).Track.cellStability > 0.5
        ratemaps.cellN{c}.session(s).Track.stableCellOn = 1;
    else
        ratemaps.cellN{c}.session(s).Track.stableCellOn = 0;
    end 
    cellStability = [cellStability, ratemaps.cellN{c}.session(s).Track.cellStability];
end 

ratemaps.vr(s).Track.PCstability = mean(cellStability, 'omitnan'); 

if ratemaps.vr(s).Track.PCstability < 0.5 
    return;
end 

disp('compute single cell correlations')
% compute spatial single cell correlations within (odd-even trials) and
% across rooms (all trials, and odd-even trials as well)
% bins per room
for r = 1:length(in.room)
    in.roomBins(r).bins = in.room(r).bins;
end 

% computing over odd vs even trials
counter = 1;
if size(ratemapsEven,1) == 1 && size(ratemapsEven,2)==sum(isnan(ratemapsEven))
    return, 
end 
for c = 1:size(ratemapsEven,1)
    cell = cellIDs(c); 
    for row = 1:2:7
        for col = 1:2:7
            if sum(ratemaps.cellN{cell}.session(s).Track.fieldPositions(1,in.roomBins(row).bins))==0 && ...
                    sum(ratemaps.cellN{cell}.session(s).Track.fieldPositions(1,in.roomBins(col).bins))==0
                continue, 
            end
            scCor.cells{counter}.session(s).rooms.oddEven(row,col) = spatialCorrelation(ratemapsOdd(c,in.roomBins(row).bins), ratemapsEven(c,in.roomBins(col).bins),0);
        end 
    end 
    for row = 2:2:7
        for col = 2:2:7
            if sum(ratemaps.cellN{cell}.session(s).Track.fieldPositions(1,in.roomBins(row).bins))==0 && ...
                    sum(ratemaps.cellN{cell}.session(s).Track.fieldPositions(1,in.roomBins(col).bins))==0
                continue, 
            end
            scCor.cells{counter}.session(s).rooms.oddEven(row,col) = spatialCorrelation(ratemapsOdd(c,in.roomBins(row).bins), ratemapsEven(c,in.roomBins(col).bins),0);
        end 
    end 
    counter = counter + 1;
end 

if isempty(scCor) 
    return, 
end 
for row = 1:7
    for col = 1:7
        vector = NaN(length(scCor.cells),1); 
        for c = 1:length(scCor.cells)
            if isempty(scCor.cells{c})
                continue, 
            end 
            if size(scCor.cells{c}.session(s).rooms.oddEven,1) < row || size(scCor.cells{c}.session(s).rooms.oddEven,2) < col
                continue, 
            end 
            if s == 3 && length(scCor.cells{c}.session)==1
                continue, 
            end 
            vector(c) = scCor.cells{c}.session(s).rooms.oddEven(row,col);
        end 
        scCor.session(s).rooms.oddEven(row,col) = mean(vector, 'omitnan'); 
   end 
end 

% now do randomized single cell correlations 
rng('shuffle')
cellSeq = randperm(size(ratemapsEven,1));
counter = 1;
for c = 1:size(ratemapsEven,1)
    random = cellSeq(c);
    randomCell = cellIDs(random); 
    cell = cellIDs(c); 
    for row = 1:2:7
        for col = 1:2:7
            if sum(ratemaps.cellN{cell}.session(s).Track.fieldPositions(1,in.roomBins(row).bins))==0 && ...
                    sum(ratemaps.cellN{randomCell}.session(s).Track.fieldPositions(1,in.roomBins(col).bins))~=0
                continue, 
            end
            scCor.randomCells{counter}.session(s).rooms.oddEven(row,col) = mean(spatialCorrelation(ratemapsOdd(c,in.roomBins(row).bins), ratemapsEven(random,in.roomBins(col).bins),0));
            
        end 
    end 
    for row = 2:2:7
        for col = 2:2:7
            if sum(ratemaps.cellN{cell}.session(s).Track.fieldPositions(1,in.roomBins(row).bins))==0 && ...
                    sum(ratemaps.cellN{randomCell}.session(s).Track.fieldPositions(1,in.roomBins(col).bins))~=0
                continue, 
            end
            scCor.randomCells{counter}.session(s).rooms.oddEven(row,col) = mean(spatialCorrelation(ratemapsOdd(c,in.roomBins(row).bins), ratemapsEven(random,in.roomBins(col).bins),0));
            
        end 
    end 
    counter = counter + 1;
end 

for row = 1:7
    for col = 1:7
       vector = NaN(length(scCor.randomCells),1); 
        for c = 1:length(scCor.randomCells)
            if isempty(scCor.randomCells{c})
                continue, 
            end 
            if size(scCor.randomCells{c}.session(s).rooms.oddEven,1) < row || size(scCor.randomCells{c}.session(s).rooms.oddEven,2) < col
                continue, 
            end 
            if s== 3 && length(scCor.randomCells{c}.session) == 1
                continue, 
            end 
            if (size(scCor.randomCells{c}.session(s).rooms.oddEven,1)<7 || size(scCor.randomCells{c}.session(s).rooms.oddEven,2)<7)
                continue, 
            end 
            vector(c) = scCor.randomCells{c}.session(s).rooms.oddEven(row,col);
        end 
        scCor.session(s).randomRooms.oddEven(row,col) = mean(vector, 'omitnan'); 
   end 
end 

% plot heatmaps of all the correlations of this session
figure()
roomLabels = {'Room A', 'Cor 1', 'Room B1', 'Cor 2', 'Room B2', 'Cor 3', 'Room C'};
h = heatmap(roomLabels,roomLabels,scCor.session(s).rooms.oddEven,'Colormap',parula, 'ColorLimits', [-1 1]); 
h.FontSize = 20;
title_str1 = sprintf ('Single cell correlations across rooms in odd vs even trials');
title_str2 = sprintf ('mouse %d day %d session %d', in.mouse, in.recordingDay, s);
title({title_str1, title_str2});
set(gcf,'units','normalized','outerposition',[0 0 1 1])
picname= sprintf('singleCellCorr_OddEven_RoomsCorridors_%d_%d_%d.png', in.mouse, in.recordingDay, s);
filepic = ['D:\Sara\Figures_analysis\Single cell correlations\' picname];
filepicname= convertCharsToStrings(filepic);
saveas(gcf,filepicname);
close(gcf)


end 