clear all

%% -- Analysis parameters

% Pixel binning, applies in all three spatial dimensions
binning = 2;

% Which color channel should be used for segmentation?
segChannel = 1;

% How many pixels should the segmentation be eroded for quantification of
% nuclear intensity?
erodeSteps = 1;

% How many pixels should the segmentation be dilated away from the initial
% segmentation to start the mask for cytoplasmic intensity quantification?
dilateSpacer = 3;
% How thick should the cytoplasmic quantification shell be?
dilateMeasure = 2;

%% --- iterative thresholding parameters

% How many steps should be used for iterative thresholding?
% At least 30 is recommended, can go up to 150 or so
threshSteps = 80;

% Parameters to connect objects between different threshold value
% segmentations. Should mostly just work like this.
distCutoff = 4; % How far can a centroid jump between two threshold values
%values before being considered a different object
volRatioCutoff = 0.5; % Change in volume to not trace object across
%thresholds, lower is more lenient, 1.0 the maximum stringency
segMinVol = 100; % minimum volume for nucleus to be recognized,
% unit: cubic micrometers
segMaxVol = 40.*10.^3; % maximum volume for nucleus to be recognized,
% unit: cubic micrometers
segMinTrackLength = 3;

% Use parallel processing?
parallel_switch = true; % allowd values: true, false

% To segment out connected components, what connectivity to use? Suggestion
% is an 18-neighbor connectivity, which means diagonal connections are
% considered as neighbors in the connected component treatment.
componentConnectivity = 18;

% Plot the nearest neighbor outputs? Leave at false!
plot_nearest_neighbor = false;

ps = parallel.Settings;
ps.Pool.AutoCreate = parallel_switch;


%% --- pick directory and find .czi files

% Requests the source directory from the user in a dialog.
sourceDir = uigetdir;

% Recursively lists all .czi files in the source directory
listing = rdir(fullfile(sourceDir,'**',filesep,'*.czi'));

numSourceFiles = numel(listing); % number of source files
sourcePaths = cell(1,numSourceFiles); % paths to the source files
filenames = cell(1,numSourceFiles); % individual file names

% Each file can contain in principle one or several stacks. So we have to
% store how many stacks are in each given file.
num_stacks_in_file = zeros(1,numSourceFiles);
num_channels_in_stack = zeros(1,numSourceFiles);


%% --- analyze image stacks

nuc_count_cell = cell(1,numSourceFiles);
NN_median_cell = cell(1,numSourceFiles);
NN_distances_cell = cell(1,numSourceFiles);

ymaxProj_cell = cell(1,numSourceFiles);
xmaxProj_cell = cell(1,numSourceFiles);
zmaxProj_cell = cell(1,numSourceFiles);

stackSize_cell = cell(1,numSourceFiles);
voxelSize_cell = cell(1,numSourceFiles);

nucInt_cell = cell(1,numSourceFiles);
cytoInt_cell = cell(1,numSourceFiles);

centroid_cell = cell(1,numSourceFiles);

vol_cell = cell(1,numSourceFiles);

if parallel_switch
    fprintf('Starting to process files in parallel.\n')
    fprintf('Progress dealing with individual files not displayed.\n')
    fprintf('Showing only percentage of files fully processed.\n')
    
    parfor_progress(numSourceFiles);
    
end


rawVoxelSizeX_vec = zeros(1,numSourceFiles);
rawVoxelSizeY_vec = zeros(1,numSourceFiles);
rawVoxelSizeZ_vec = zeros(1,numSourceFiles);

for ff = 1:numSourceFiles
    
    thisPath = listing(ff).name;
    
    % --- Make a reader instance
    reader = bfGetReader(thisPath);
    
    % --- extract stack and microscope info from meta data
    omeMeta = reader.getMetadataStore();
    
    reader.close();
    
    numChannels = omeMeta.getChannelCount(0);
    numImages = omeMeta.getImageCount();
    
    voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0);
    voxelSizeX = voxelSizeX.value(ome.units.UNITS.MICROM);
    rawVoxelSizeX_vec(ff) = voxelSizeX.doubleValue();
    voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0);
    voxelSizeY = voxelSizeY.value(ome.units.UNITS.MICROM);
    rawVoxelSizeY_vec(ff) = voxelSizeY.doubleValue();
    voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0);
    voxelSizeZ = voxelSizeZ.value(ome.units.UNITS.MICROM);
    rawVoxelSizeZ_vec(ff) = voxelSizeZ.doubleValue();
    
end

errorFlagVec = cell(1,numSourceFiles);
errorMessages = cell(1,numSourceFiles);

parfor ff = 1:numSourceFiles
    
    thisPath = listing(ff).name;
    sourcePaths{ff} = thisPath;
    
    if ~parallel_switch
        
        fprintf('Accessing file %d of %d files (%s):\n',...
            ff,numSourceFiles,sourcePaths{ff})
        
    end
    
    % --- Make a reader instance
    reader = bfGetReader(thisPath);
    
    % --- extract stack and microscope info from meta data
    omeMeta = reader.getMetadataStore();
    
    reader.close();
    
    numChannels = omeMeta.getChannelCount(0);
    numImages = omeMeta.getImageCount();
    
    if ~parallel_switch
        fprintf('File contains %d stacks, with %d color channels.\n',...
            numImages,numChannels)
    end
    
    num_stacks_in_file(ff) = numImages;
    
    num_channels_in_stack(ff) = numChannels;
    
    rawStackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
    rawStackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
    rawStackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices
    
    if ~parallel_switch
        fprintf('Raw stack dimensions: %d by %d by %d voxels.\n',...
            rawStackSizeY,rawStackSizeX,rawStackSizeZ)
    end
    
    centerSliceInd = round(rawStackSizeX./2);
    
    rawVoxelSizeX = rawVoxelSizeX_vec(ff);
    rawVoxelSizeY = rawVoxelSizeY_vec(ff);
    rawVoxelSizeZ = rawVoxelSizeZ_vec(ff);
    
    rawVoxelVol = rawVoxelSizeX.*rawVoxelSizeY.*rawVoxelSizeZ;
    
    errorFlagVec{ff} = false(1,numImages);
    errorMessages{ff} = cell(1,numImages);
    
    numNuc = zeros(1,numImages);
    nucVol = cell(1,numImages);
    nucCent = cell(1,numImages);
    nucBBox = cell(1,numImages);
    nucVxlIdx = cell(1,numImages);
    nucImg = cell(1,numImages);
    
    nucInt = cell(1,numImages);
    cytoInt = cell(1,numImages);
    
    % Store maximum z projection and central x section
    maxProjCellYY = cell(1,numImages);
    maxProjCellXX = cell(1,numImages);
    maxProjCellZZ = cell(1,numImages);
    nucProjCell = cell(1,numImages);
    xNucSectionCell = cell(1,numImages);
    
    stackSizes = cell(1,numImages);
    voxelSizes = cell(1,numImages);
    
    % Nearest neighbor distance containers
    NN_distances = cell(1,numImages);
    NN_median = zeros(1,numImages);
    NN_std = zeros(1,numImages);
    
    
    for kk = 1:numImages
        
        if ~parallel_switch
            fprintf('Processing stack %d/%d of group %d/%d\n',...
                kk,numImages,ff,numSourceFiles);
        end
        
        
        try % Prevents crash when something is wrong with a specific file
           
            reader = bfGetReader(thisPath);

            reader.setSeries(kk-1);
            
            % Read in the raw stack
            
            rawStack = cell(1,numChannels);
            
            if ~parallel_switch
                fprintf('Reading in stack.\n')
                
                parfor_progress(numChannels.*rawStackSizeZ);
            end
            
            for cc = 1:numChannels
                
                rawStack{cc} = ...
                    zeros(rawStackSizeY,rawStackSizeX,rawStackSizeZ,'uint16');
                
                for ll = 1:rawStackSizeZ
                    
                    % Direct assigment
                    planeInd = reader.getIndex(ll-1,cc-1,0)+1;
                    rawStack{cc}(:,:,ll) = bfGetPlane(reader,planeInd);
                    
                    if ~parallel_switch
                        parfor_progress;
                    end
                    
                end
                
            end
            
            if ~parallel_switch
                parfor_progress(0);
            end
            reader.close();
            
            
            % --- Make binned stack
            if binning > 1
                
                if ~parallel_switch
                    fprintf('Binning stack, %d binning size.\n',binning)
                end
                
                binnedStack = cell(1,numChannels);
                
                [m,n,o]=size(rawStack{1});
                
                m = binning.*floor(m/binning);
                n = binning.*floor(n/binning);
                o = binning.*floor(o/binning);
                
                if ~parallel_switch
                    parfor_progress(numChannels.*o);
                end
                
                for cc = 1:numChannels
                    
                    % dimension 1 binning
                    
                    xyBinStack = zeros(m/binning,n/binning,o/binning);
                    
                    for ll = 1:o
                        
                        container = squeeze(rawStack{cc}(1:m,1:n,ll));
                        
                        container = sum(reshape(container,binning,[]),1);
                        container = reshape(container,m/binning,[]);
                        container = permute(container,[2,1]);
                        
                        container = sum(reshape(container,binning,[]),1);
                        container = reshape(container,n/binning,[]);
                        
                        xyBinStack(:,:,ll) = permute(container,[2,1]);
                        
                        if ~parallel_switch
                            parfor_progress;
                        end
                        
                    end
                    
                    xyBinStack = permute(xyBinStack,[3,2,1]);
                    xyBinStack = sum(reshape(xyBinStack,binning,[],m/binning),1);
                    xyBinStack = reshape(xyBinStack,o/binning,[],m/binning);
                    xyBinStack = permute(xyBinStack,[3,2,1]);
                    
                    binnedStack{cc} = xyBinStack./binning.^3;
                    
                end
                
                if ~parallel_switch
                    parfor_progress(0);
                end
                
            else
                
                binnedStack = rawStack;
                
            end
            
            rawStack = 0;
            
            [stackSizeY,stackSizeX,stackSizeZ] = size(binnedStack{segChannel});
            voxelSizeY = rawVoxelSizeY.*binning;
            voxelSizeX = rawVoxelSizeY.*binning;
            voxelSizeZ = rawVoxelSizeZ.*binning;
            voxelVol = voxelSizeY.*voxelSizeX.*voxelSizeZ;
            
            maxProjCellYY{kk} = squeeze(max(binnedStack{segChannel},[],1));
            maxProjCellXX{kk} = squeeze(max(binnedStack{segChannel},[],2));
            maxProjCellZZ{kk} = squeeze(max(binnedStack{segChannel},[],3));
            
            voxelSizes{kk} = [voxelSizeY,voxelSizeX,voxelSizeZ];
            stackSizes{kk} = [stackSizeY,stackSizeX,stackSizeZ];
            
            volCell = cell(1,threshSteps);
            centrCell = cell(1,threshSteps);
            vxlIdxCell = cell(1,threshSteps);
            coordsVxlIdxCell = cell(1,threshSteps);
            imagesCell = cell(1,threshSteps);
            
            numRegVec = zeros(1,threshSteps);
            
            featureVecCell = cell(1,threshSteps);
            
            if ~parallel_switch
                fprintf('Nuclei segmentation with iterative thresholds.\n')
            end
            
            % --- get threshold values from percentile steps
            
            jointInts = [maxProjCellYY{kk}(:);...
                maxProjCellXX{kk}(:);...
                maxProjCellZZ{kk}(:)];
            
            [otsuThresh,~] = otsuLimit(jointInts);
            
            threshVals = ...
                linspace(otsuThresh.*0.5,max(jointInts),threshSteps+1);
            
            %         scalingExponent = 3;
            %         threshVals = linspace(...
            %             (otsuThresh.*0.5).^(-scalingExponent), ...
            %             max(jointInts).^(-scalingExponent),threshSteps+1)...
            %             .^scalingExponent;
            
            threshVals = threshVals(1:end-1);
            
            if ~parallel_switch
                parfor_progress(threshSteps);
            end
            
            for tt = 1:threshSteps
                
                currThr = threshVals(tt);
                
                regions = bwconncomp(binnedStack{segChannel}>=currThr,...
                    componentConnectivity);
                regionProps = regionprops(regions,'Area');
                regionVols = [regionProps(:).Area].*voxelVol;
                
                % Restrict to objects greater than minimum volume
                validVolInds = regionVols>segMinVol & regionVols<segMaxVol;
                nucleiRegions = struct;
                nucleiRegions.Connectivity = regions.Connectivity;
                nucleiRegions.ImageSize = regions.ImageSize;
                nucleiRegions.NumObjects = sum(validVolInds);
                nucleiRegions.PixelIdxList = regions.PixelIdxList(validVolInds);
                
                regionProps = ...
                    regionprops(nucleiRegions,'Area','Centroid',...
                    'PixelIdxList');
                volumes = [regionProps(:).Area].*voxelVol;
                centroids = {regionProps(:).Centroid};
                vxlIdx = {regionProps(:).PixelIdxList};
                
                yCoords = cellfun(@(elmt)elmt(1).*voxelSizeY,centroids);
                xCoords = cellfun(@(elmt)elmt(2).*voxelSizeX,centroids);
                zCoords = cellfun(@(elmt)elmt(3).*voxelSizeZ,centroids);
                
                featureVecCell{tt} = [volumes;xCoords;yCoords;zCoords].';
                
                numRegVec(tt) = numel(volumes);
                volCell{tt} = volumes;
                centrCell{tt} = centroids;
                vxlIdxCell{tt} = vxlIdx;
                
                if ~parallel_switch
                    parfor_progress;
                end
                
            end
            
            if ~parallel_switch
                parfor_progress(0);
            end
            
            % --- Linking of objects
            
            if ~parallel_switch
                fprintf('Selection of individual nuclei thresholds.\n')
                
                parfor_progress(2.*threshSteps);
            end
            rvsLinkInds = cell(1,threshSteps-1);
            physicalDist = cell(1,threshSteps-1);
            volRatio = cell(1,threshSteps-1);
            
            for tt = 1:(threshSteps-1)
                
                if ~isempty(featureVecCell{tt+1}) ...
                        && ~isempty(featureVecCell{tt})
                    
                    distMatr = ...
                        pdist2(featureVecCell{tt+1}(:,2:end),...
                        featureVecCell{tt}(:,2:end),'euclidean');
                    
                    [rvsLinkDist,rvsLinkInds{tt}] = min(distMatr,[],2);
                    
                    physicalDist{tt} = sqrt(sum((...
                        featureVecCell{tt}(rvsLinkInds{tt},2:end) ...
                        - featureVecCell{tt+1}(:,2:end)).^2,2));
                    
                    volRatio{tt} = ...
                        featureVecCell{tt+1}(:,1)...
                        ./featureVecCell{tt}(rvsLinkInds{tt},1);
                    volRatio{tt}(volRatio{tt}>1) = 1./volRatio{tt}(volRatio{tt}>1);
                    
                else
                    rvsLinkInds{tt} = [];
                    physicalDist{tt} = [];
                    volRatio{tt} = [];
                end
                
                if ~parallel_switch
                    parfor_progress;
                end
                
            end
            
            objects = {};
            if numRegVec(end)>0
                for rr = 1:numRegVec(end)
                    
                    thisRegion = struct;
                    
                    thisRegion.origObjIndTrack = rr;
                    thisRegion.volTrack = volCell{end}(rr);
                    thisRegion.threshTrack = threshSteps;
                    thisRegion.centrTrack = centrCell{end}(rr);
                    thisRegion.truncated = false;
                    
                    objects = [objects,{thisRegion}];
                    
                end
            end
            
            
            for tt = (threshSteps-1):-1:1
                
                nonTruncInds = ...
                    find(cellfun(@(elmt)~elmt.truncated,objects));
                
                usedInds = [];
                
                if ~isempty(rvsLinkInds{tt})
                    
                    % make all objects truncated, and only bring those back that
                    % get connected
                    
                    for oo = 1:numel(objects)
                        objects{oo}.truncated = true;
                    end
                    
                    useFlags = ...
                        (volRatio{tt}>=volRatioCutoff ...
                        & physicalDist{tt}<=distCutoff);
                    
                    %                 disp('+++')
                    %                 disp(numRegVec(tt+1))
                    %                 disp(numel(nonTruncInds))
                    
                    for rr = 1:numRegVec(tt+1)
                        
                        objInd = nonTruncInds(rr);
                        origObjInd = objects{objInd}.origObjIndTrack(end);
                        rvsInd = rvsLinkInds{tt}(origObjInd);
                        
                        if useFlags(origObjInd)
                            
                            thisRegion = struct;
                            
                            objects{objInd}.origObjIndTrack = ...
                                [objects{objInd}.origObjIndTrack,rvsInd];
                            objects{objInd}.volTrack = ...
                                [objects{objInd}.volTrack,volCell{tt}(rvsInd)];
                            objects{objInd}.threshTrack = ...
                                [objects{objInd}.threshTrack,tt];
                            objects{objInd}.centrTrack = ...
                                [objects{objInd}.centrTrack,centrCell{tt}(rvsInd)];
                            objects{objInd}.truncated = false;
                            
                            usedInds = [usedInds,rvsInd];
                            
                        end
                        
                    end
                    
                end
                
                %             disp(numel(usedInds))
                notUsedInds = setdiff(1:numRegVec(tt),usedInds);
                %             disp(numel(notUsedInds))
                
                
                
                for rr = 1:numel(notUsedInds)
                    
                    regInd = notUsedInds(rr);
                    
                    thisRegion = struct;
                    
                    thisRegion.origObjIndTrack = regInd;
                    thisRegion.volTrack = volCell{tt}(regInd);
                    thisRegion.threshTrack = tt;
                    thisRegion.centrTrack = centrCell{tt}(regInd);
                    thisRegion.truncated = false;
                    
                    objects = [objects,{thisRegion}];
                    
                end
                
                if ~parallel_switch
                    parfor_progress;
                end
                
            end
            
            
            % -- Find individual objects' threshold by minimum volume change
            
            keepFlags = ...
                cellfun(@(elmt)numel(elmt.threshTrack)>segMinTrackLength,objects);
            keepObjects = objects(keepFlags);
            
            for oo = 1:numel(keepObjects)
                
                %     volDiff = diff(keepObjects{oo}.volTrack) ...
                %         ./keepObjects{oo}.volTrack(1:end-1);
                %
                %
                centrs = keepObjects{oo}.centrTrack;
                scalVec = [voxelSizeY,voxelSizeX,voxelSizeZ];
                centrDiff = zeros(1,numel(centrs)-1);
                for pp = 1:numel(centrs)-1
                    centrDiff(pp) = sqrt( ...
                        sum((centrs{pp}.*scalVec-centrs{pp+1}.*scalVec).^2));
                end
                centrDiff = centrDiff./keepObjects{oo}.volTrack(1:end-1);
                minDiff = min(centrDiff);
                minDiffInd = find(centrDiff==minDiff,1,'last');
                
                %             minDiffInd = numel(keepObjects{oo}.volTrack)-1;
                keepObjects{oo}.threshInd = keepObjects{oo}.threshTrack(minDiffInd);
                keepObjects{oo}.centr = keepObjects{oo}.centrTrack{minDiffInd};
                keepObjects{oo}.vol = keepObjects{oo}.volTrack(minDiffInd);
                keepObjects{oo}.origObjInd = ...
                    keepObjects{oo}.origObjIndTrack(minDiffInd);
                keepObjects{oo}.vxlIdx = ...
                    vxlIdxCell{keepObjects{oo}.threshInd}{keepObjects{oo}.origObjInd};
                
            end
            
            if ~parallel_switch
                parfor_progress;
            end
            
            % -- rejection of overlapping objects
            
            numObjs = numel(keepObjects);
            keepFlags = true(1,numObjs);
            
            yCoords = cellfun(@(elmt)elmt.centr(1).*voxelSizeY,keepObjects);
            xCoords = cellfun(@(elmt)elmt.centr(2).*voxelSizeX,keepObjects);
            zCoords = cellfun(@(elmt)elmt.centr(3).*voxelSizeZ,keepObjects);
            
            coordVector = [xCoords.',yCoords.',zCoords.'];
            objectDist = pdist2(coordVector,coordVector);
            
            for rr = (numRegVec(end)+1):numObjs
                
                %             numObjs - rr
                
                checkInds = 1:(rr-1);
                [dists,sortInds] = sort(objectDist(rr,checkInds),'ascend');
                
                sortInds = sortInds(dists<50); % limit to a sphere of 50 microns radius
                
                checkInds = checkInds(sortInds);
                
                for qq = checkInds;
                    
                    if ~isempty(intersect(keepObjects{rr}.vxlIdx,keepObjects{qq}.vxlIdx))
                        
                        keepFlags(rr) = false;
                        break;
                        
                    end
                    
                end
                
            end
            
            validObjs = keepObjects(keepFlags);
            
            % -- Make binary image of valid objects
            % Add all objects to binary 3D image
            segArray = false(stackSizeY,stackSizeX,stackSizeZ);
            numValObj = numel(validObjs);
            for oo = 1:numValObj
                segArray(validObjs{oo}.vxlIdx) = true;
            end
            
            % Store maximum z projection and central x section
            nucProjCell{kk} = max(segArray,[],3);
            xNucSectionCell{kk} = squeeze(segArray(:,centerSliceInd,:));
            
            if ~parallel_switch
                parfor_progress;
                parfor_progress(0);
            end
            
            % -- Quantification based on extracted regions
            
            if ~parallel_switch
                fprintf('Quantification of nuclear and cytoplasmic intensities.\n')
            end
            
            % Construct region structure for further intensity feature
            % extraction
            nucleiRegions = struct;
            nucleiRegions.Connectivity = componentConnectivity;
            nucleiRegions.ImageSize = [stackSizeY,stackSizeX,stackSizeZ];
            nucleiRegions.NumObjects = numel(validObjs);
            nucleiRegions.PixelIdxList = ...
                cellfun(@(elmt)elmt.vxlIdx,validObjs,'UniformOutput',false);
            
            nucleiPropsBinary = regionprops(nucleiRegions,...
                'Area','Image','BoundingBox');
            
            nucleiVol = [nucleiPropsBinary(:).Area].*voxelVol;
            nucleiImage = {nucleiPropsBinary(:).Image};
            nucleiBoundingBox = {nucleiPropsBinary(:).BoundingBox};
            nucleiVxlIdx = nucleiRegions.PixelIdxList;
            
            nucleiPropsSegmentation = ...
                regionprops(nucleiRegions,binnedStack{segChannel},...
                'WeightedCentroid');
            
            nucleiCentroid = ...
                {nucleiPropsSegmentation(:).WeightedCentroid};
            
            numNuc(kk) = numel(nucleiVol);
            nucVol{kk} = nucleiVol;
            nucCent{kk} = nucleiCentroid;
            for ll = 1:numNuc(kk)
                nucCent{kk}{ll} = ...
                    nucCent{kk}{ll}.*[voxelSizeY,voxelSizeX,voxelSizeZ];
            end
            nucBBox{kk} = nucleiBoundingBox;
            nucVxlIdx{kk} = nucleiVxlIdx;
            nucImg{kk} = nucleiImage;
            
            nucInt{kk} = cell(1,numChannels);
            cytoInt{kk} = cell(1,numChannels);
            
            for cc = 1:numChannels
                
                nucInt{kk}{cc} = zeros(1,numNuc(kk));
                cytoInt{kk}{cc} = zeros(1,numNuc(kk));
                
            end
            
            
            % Make masks to determine nuclear and cytoplasmic intensities
            
            totalDil = dilateSpacer+dilateMeasure;
            
            if ~parallel_switch
                parfor_progress(numNuc(kk));
            end
            
            for nn = 1:numNuc(kk)
                
                % Determine small bounding box (without dilate)
                smallMinY = nucBBox{kk}{nn}(2)+0.5;
                smallMinX = nucBBox{kk}{nn}(1)+0.5;
                smallMinZ = nucBBox{kk}{nn}(3)+0.5;
                smallMaxY = nucBBox{kk}{nn}(2)+nucBBox{kk}{nn}(5)-0.5;
                smallMaxX = nucBBox{kk}{nn}(1)+nucBBox{kk}{nn}(4)-0.5;
                smallMaxZ = nucBBox{kk}{nn}(3)+nucBBox{kk}{nn}(6)-0.5;
                
                
                % Determine extended bounding box (after dilate)
                fullExtMinY = smallMinY-totalDil;
                fullExtMinX = smallMinX-totalDil;
                fullExtMinZ = smallMinZ-totalDil;
                fullExtMaxY = smallMaxY+totalDil;
                fullExtMaxX = smallMaxX+totalDil;
                fullExtMaxZ = smallMaxZ+totalDil;
                
                % Limit extended bounding box to within image limits
                extMinY = max(1,fullExtMinY);
                yLoDiff = extMinY - fullExtMinY;
                extMinX = max(1,fullExtMinX);
                xLoDiff = extMinX - fullExtMinX;
                extMinZ = max(1,fullExtMinZ);
                zLoDiff = extMinZ - fullExtMinZ;
                extMaxY = min(stackSizeY,fullExtMaxY);
                yHiDiff = fullExtMaxY - extMaxY;
                extMaxX = min(stackSizeX,fullExtMaxX);
                xHiDiff = fullExtMaxX - extMaxX;
                extMaxZ = min(stackSizeZ,fullExtMaxZ);
                zHiDiff = fullExtMaxZ - extMaxZ;
                
                % Extended bounding box size
                extSizeY = extMaxY - extMinY + 1;
                extSizeX = extMaxX - extMinX + 1;
                extSizeZ = extMaxZ - extMinZ + 1;
                
                % Nucleus mask, will be eroded into nucleus
                nuclMask = false(extSizeY,extSizeX,extSizeZ);
                nuclMask((1+totalDil-yLoDiff):(end-totalDil+yHiDiff),...
                    (1+totalDil-xLoDiff):(end-totalDil+xHiDiff),...
                    (1+totalDil-zLoDiff):(end-totalDil+zHiDiff))...
                    = nucImg{kk}{nn};
                
                % Inclusion mask for cytoplasm shell, growing from nucleus
                % under investigation only
                inclMask = nuclMask;
                
                % Exclusion mask for cytoplasm shell, growing from all
                % segmented objects
                exclMask = segArray(extMinY:extMaxY,...
                    extMinX:extMaxX,...
                    extMinZ:extMaxZ);
                
                dilateNeighborhood = zeros(3,3,3);
                dilateNeighborhood(2,2,2) = 1;
                dilateNeighborhood(1,2,2) = 1;
                dilateNeighborhood(3,2,2) = 1;
                dilateNeighborhood(2,1,2) = 1;
                dilateNeighborhood(2,3,2) = 1;
                dilateNeighborhood(2,2,1) = 1;
                dilateNeighborhood(2,2,3) = 1;
                
                for ee = 1:erodeSteps
                    nuclMask = imerode(nuclMask,dilateNeighborhood);
                end
                
                for dd = 1:dilateSpacer
                    inclMask = imdilate(inclMask,dilateNeighborhood);
                    exclMask = imdilate(exclMask,dilateNeighborhood);
                end
                
                for dd = 1:dilateMeasure
                    inclMask = imdilate(inclMask,dilateNeighborhood);
                end
                
                measureMask = inclMask & ~exclMask;
                
                for cc = 1:numChannels
                    
                    % -- Determine nuclear to cytoplasmic ratios
                    
                    cutoutImage = binnedStack{cc}(extMinY:extMaxY,...
                        extMinX:extMaxX,...
                        extMinZ:extMaxZ);
                    
                    nucInt{kk}{cc}(nn) = mean(cutoutImage(nuclMask));
                    cytoInt{kk}{cc}(nn) = mean(cutoutImage(measureMask));
                    
                end
                
                if ~parallel_switch
                    parfor_progress;
                end
                
            end
            
            if ~parallel_switch
                parfor_progress(0);
            end
            
            % --- Nearest neighbor distances
            
            if ~parallel_switch
                fprintf('Calculating nearest neighbor distances...')
            end
            
            centroids = zeros(numNuc(kk),3);
            
            for nn = 1:numNuc(kk)
                centroids(nn,:) = nucCent{kk}{nn};
            end
            
            % calculate median pairwise distance to nearest neighbor
            distMatrix = squareform(pdist(centroids));
            distMatrix(distMatrix == 0) = Inf;
            [distances,NNinds] = min(distMatrix);
            
            NN_distances{kk} = distances;
            NN_median(kk) = median(distances);
            NN_std(kk) = std(distances);
            
            if plot_nearest_neighbor && ~parallel_switch
                
                % --- plotting of nearest neighbor distances
                
                figure(2)
                
                plot(NN_median(kk).*[1,1],[0,1],'r--')
                
                hold on
                
                plot(sort(distances),linspace(0,1,numel(distances)),'k-o')
                
                hold off
                
                xlabel('Individual nucleus NN distance [\mum]')
                ylabel('Cumulative probability')
                
            end
            
            if ~parallel_switch
                fprintf('done.\n')
            end
            
            % Upon successful completion, set error flag to false, no error!
            errorFlagVec{ff}(kk) = false;
<<<<<<< HEAD
                        
=======

>>>>>>> origin/master
        catch message
            
            errorFlagVec{ff}(kk) = true;
            errorMessages{ff}{kk} = message;
            
            fprintf('Image segmentation error in step %d.\n',kk)
            
        end
        
    end
    
    nuc_count_cell{ff} = numNuc;
    NN_distances_cell{ff} = NN_distances;
    NN_median_cell{ff} = NN_median;
    
    ymaxProj_cell{ff} = maxProjCellYY;
    xmaxProj_cell{ff} = maxProjCellXX;
    zmaxProj_cell{ff} = maxProjCellZZ;
    
    stackSize_cell{ff} = stackSizes;
    voxelSize_cell{ff} = voxelSizes;
    
    nucInt_cell{ff} = nucInt;
    cytoInt_cell{ff} = cytoInt;
    
    centroid_cell{ff} = nucCent;
    
    vol_cell{ff} = nucVol;
    
    if parallel_switch
        parfor_progress;
    end
    
end

parfor_progress(0);

%% --- save individual embryo results for further processing

fprintf('Saving individual embryo results...\n')

parfor_progress(numSourceFiles);

for ff = 1:numSourceFiles
    
    parfor_progress;
    
    saveStruct = struct;
    
    numEmbryos = numel(nuc_count_cell{ff});
    
    for ee = 1:numEmbryos
        
        if ~errorFlagVec{ff}(ee)
            
            saveStruct.nuc_count = nuc_count_cell{ff}(ee);
            saveStruct.NN_median = NN_median_cell{ff}(ee);
            saveStruct.NN_distances = NN_distances_cell{ff}{ee};
            
            saveStruct.ymaxProj = ymaxProj_cell{ff}{ee};
            saveStruct.xmaxProj = xmaxProj_cell{ff}{ee};
            saveStruct.zmaxProj = zmaxProj_cell{ff}{ee};
            
            saveStruct.stackSize = stackSize_cell{ff}{ee};
            saveStruct.voxelSize = voxelSize_cell{ff}{ee};
            
            saveStruct.nucInt = nucInt_cell{ff}{ee};
            saveStruct.cytoInt = cytoInt_cell{ff}{ee};
            
            saveStruct.centroid = centroid_cell{ff}{ee};
            saveStruct.volume = vol_cell{ff}{ee};
            
            [pathOnly,nameOnly,extensionOnly] = fileparts(sourcePaths{ff});
            
            mkdir(fullfile(pathOnly,filesep,'Results'))
            
            saveFileName = fullfile(...
                pathOnly,filesep,'Results',filesep,...
                sprintf('analyzed_%d_%d.mat',ff,ee));
                        
            save(saveFileName,'-struct','saveStruct');
            
        end
        
    end
    
end

parfor_progress(0);

fprintf('Saving complete.\n')

%% --- plot overview of nuclei counts and nearest neighbor distances

joint_dist_medians = [NN_median_cell{:}];
joint_nuc_counts = [nuc_count_cell{:}];

joint_dist_distrs = [NN_distances_cell{:}];

upper_lim = cellfun(@(distr) prctile(distr,60),joint_dist_distrs);
lower_lim = cellfun(@(distr) prctile(distr,40),joint_dist_distrs);


errorbar(joint_nuc_counts,joint_dist_medians,...
    lower_lim-joint_dist_medians,upper_lim-joint_dist_medians,...
    'ko')

ylabel('Nearest neighbor distance [\mum]')
xlabel('Nuclei count')
