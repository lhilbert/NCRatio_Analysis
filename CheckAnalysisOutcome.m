clear all


%% --- analysis parameters

% Put in channel number to use as marker for interphase
% (e.g. Pol II or PCNA)
interphaseMarkerChannel = 1;

% Put in channel number from which nuc/cyto ratio intensities should be put
% out for intensity investigation
intChannel = 2;

% Minimum nuc/cyto intensity ratio of the interphase marker to include
% nucleus in neigbor distance analysis
minIntRatio = 2;

% Minimum number of nuclei to make valid conclusion about neighbor distance
minNucCount = 15;

% Do you want to see the output of every single file? assign true or false
showResults = true;

%% --- pick directory and find all files

sourceDir = uigetdir;

listing = dir(fullfile(sourceDir,'*.mat'));
numResultFiles = numel(listing);
sourcePaths = cell(1,numResultFiles);
filenames = cell(1,numResultFiles);


%% --- run through files and show maximum projections with marked nuclei

% Write file opening
writeFile = fopen(fullfile(sourceDir,'stagingTable.dat'),'w+');

fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('Source directory: %s\n',sourceDir)
fprintf('Embryo file, Median distance [um], Interphase nuclei\n',...
    sourceDir);

fprintf(writeFile,'Source directory: %s\n',sourceDir);
fprintf(writeFile,'Embryo file, Median distance [um], Interphase nuclei\n',...
    sourceDir);


medianDistanceVec = zeros(1,numResultFiles);
numValidNuc = zeros(1,numResultFiles);

for ff = 1:numResultFiles
    
    try
    
    figure(1)
    
    clf
    
    % load result file
    
    loadedStruct = load([sourceDir,filesep,listing(ff).name]);
    
    xxCoords = cellfun(@(elmt)elmt(1),loadedStruct.centroid);
    yyCoords = cellfun(@(elmt)elmt(2),loadedStruct.centroid);
    zzCoords = cellfun(@(elmt)elmt(3),loadedStruct.centroid);
    
    
    % -- Get nearest neighbor distance only for sufficiently bright nuclei
    
    % Indices of nuclei with high enough intensity
    
    if loadedStruct.nuc_count>0
        
        nucInts = loadedStruct.nucInt{interphaseMarkerChannel};
        cytoInts = loadedStruct.cytoInt{interphaseMarkerChannel};
        
        intRatios = nucInts./cytoInts;
        
        validNucInds = find(intRatios>=minIntRatio);
        
        numValidNuc(ff) = numel(validNucInds);
        
        selectedNucleiCentroids = ...
            loadedStruct.centroid(validNucInds);
        
    else
        
        numValidNuc(ff) = 0;
        selectedNucleiCentroids = {};
        
    end
    
    
    
    
    % --- Nearest neighbor distances
    
    if numValidNuc(ff) <= minNucCount
        
        distances = [];
        median_val = NaN;
        var_val = NaN;
        
    else

        xxCoords = cellfun(@(elmt)elmt(1),selectedNucleiCentroids);
        yyCoords = cellfun(@(elmt)elmt(2),selectedNucleiCentroids);
        zzCoords = cellfun(@(elmt)elmt(3),selectedNucleiCentroids);

        coord_matrix = [yyCoords;xxCoords;zzCoords].';
        
        % calculate median pairwise distance to second-nearest neighbor
        distMatrix = squareform(pdist(coord_matrix));
        distMatrix(distMatrix == 0) = Inf;
        distMatrix = sort(distMatrix);
        distances = squeeze(distMatrix(2,:));
        
        median_val = median(distances);
        var_val = var(distances);
        
    end
        
    if showResults
        
        subplot(2,3,1)
        
        imagesc([0,loadedStruct.stackSize(1).*loadedStruct.voxelSize(1)],...
            [0,loadedStruct.stackSize(2).*loadedStruct.voxelSize(2)],...
            loadedStruct.zmaxProj)
        
        hold on
        
        plot(xxCoords,yyCoords,'r+')
        
        xlabel('x [\mum]')
        ylabel('y [\mum]')
        
        colormap(gray)
        
        axis equal
        axis tight
        
        titleHandle = title(sprintf('File name: %s',listing(ff).name));
        set(titleHandle,'interpreter','none')
        
        
        subplot(2,3,2)
        
        imagesc([0,loadedStruct.stackSize(2).*loadedStruct.voxelSize(2)],...
            [0,loadedStruct.stackSize(3).*loadedStruct.voxelSize(3)],...
            loadedStruct.xmaxProj.')
        
        hold on
        
        plot(yyCoords,zzCoords,'r+')
        
        xlabel('y [\mum]')
        ylabel('z [\mum]')
        
        colormap(gray)
        
        axis equal
        axis tight
        
        set(gca,'YDir','normal')

        titleHandle = title(sprintf('Number of colors: %d',...
            numel(loadedStruct.cytoInt)));
        set(titleHandle,'interpreter','none')

        
        subplot(3,3,3)
        
        imagesc([0,loadedStruct.stackSize(1).*loadedStruct.voxelSize(1)],...
            [0,loadedStruct.stackSize(3).*loadedStruct.voxelSize(3)],...
            loadedStruct.ymaxProj.')
        
        hold on
        
        plot(xxCoords,zzCoords,'r+')
        
        xlabel('y [\mum]')
        ylabel('z [\mum]')
        
        colormap(gray)
        
        axis equal
        axis tight
        
        set(gca,'YDir','normal')
        
        
        % --- plot nuclear intensity histogram
        
        subplot(2,3,4)
        plot(loadedStruct.cytoInt{intChannel},...
            loadedStruct.nucInt{intChannel},'ko')
        hold on
        minMaxSupport = [0,max(loadedStruct.nucInt{intChannel})];
        plot(minMaxSupport,minMaxSupport,'k--')
        xlabel('Cytoplasmic intensity')
        ylabel('Nuclear intensity')
        
        subplot(2,3,5)
        hist(loadedStruct.nucInt{intChannel} ...
            ./loadedStruct.cytoInt{intChannel},20)
        xlabel('Intensity ratio Nuc/Cyto')
        ylabel('Count')
        
        
        % --- plot neighbor distance distribution
        
        
        subplot(2,3,6)
        hist(distances,20)
        xlabel('Second neighbor distance [\mum]')
        ylabel('Count')
        
        waitforbuttonpress
        
    end

    % --- print out median neighbor distance for different files
    
    fprintf('%s,%6.6f,%d\n',listing(ff).name,median_val,numValidNuc(ff))

    fprintf(writeFile,'%s,%6.6f,%d\n',listing(ff).name,median_val,numValidNuc(ff));

    catch ME
        
        disp('Failure in file.')
        
    end
    
    
end

fclose(writeFile);