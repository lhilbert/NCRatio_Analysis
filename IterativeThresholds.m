clear all

[filepath,folderpath] = uigetfile('*.*');


%% -- Analysis parameters

componentConnectivity = 18;

binning = 2;

% Assign channels for different purposes
segChannel = 1; % Segmentation of nuclei (DAPI)

dilateSpacer = 3;
dilateMeasure = 4;

% --- iterative thresholding parameters

threshSteps = 100;
maxThresh = 40000;
minThresh = 500;

threshVals = linspace(minThresh,maxThresh,threshSteps);

distCutoff = 4; % How far can a centroid jump between two threshold values
%values before being considered a different object
volRatioCutoff = 0.5; % Change in volume to not trace object across
%thresholds, lower is more lenient, 1.0 the maximum stringency
segMinVol = 50; % minimum volume for nucleus to be recognized,
% unit: cubic micrometers
segMaxVol = 30.*10.^3; % maximum volume for nucleus to be recognized,
% unit: cubic micrometers
segMinTrackLength = 3;



plot_nearest_neighbor = true;



% --- Make a reader instance
reader = bfGetReader([folderpath,filesep,filepath]);

% --- extract stack and microscope info from meta data
omeMeta = reader.getMetadataStore();

reader.close();

numChannels = omeMeta.getChannelCount(0);
numImages = omeMeta.getImageCount();

rawStackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
rawStackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
rawStackSizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

centerSliceInd = round(rawStackSizeX./2);

voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0);
voxelSizeX = voxelSizeX.value(ome.units.UNITS.MICROM);
rawVoxelSizeX = voxelSizeX.doubleValue();
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0);
voxelSizeY = voxelSizeY.value(ome.units.UNITS.MICROM);
rawVoxelSizeY = voxelSizeY.doubleValue();
voxelSizeZ = omeMeta.getPixelsPhysicalSizeZ(0);
voxelSizeZ = voxelSizeZ.value(ome.units.UNITS.MICROM);
rawVoxelSizeZ = voxelSizeZ.doubleValue();

rawVoxelVol = rawVoxelSizeX.*rawVoxelSizeY.*rawVoxelSizeZ;


errorFlagVec = false(1,numImages);

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


% Nearest neighbor distance containers
NN_distances{kk} = distances;
NN_median = zeros(1,kk);
NN_std = zeros(1,kk);


for kk = 1:numImages
    
    reader = bfGetReader([folderpath,filesep,filepath]);
    
    %     try
    
    reader.setSeries(kk-1);
    
    % Read in the raw stack
    
    rawStack = cell(1,numChannels);
    
    for cc = 1:numChannels
        
        rawStack{cc} = ...
            zeros(rawStackSizeY,rawStackSizeX,rawStackSizeZ,'uint16');
        
        for ll = 1:rawStackSizeZ
            
            ll
            
            % Direct assigment
            planeInd = reader.getIndex(ll-1,cc-1,0)+1;
            rawStack{cc}(:,:,ll) = bfGetPlane(reader,planeInd);
            
        end
        
    end
    
    reader.close();
    
    % --- Make binned stack
    if binning > 1
        
        binnedStack = cell(1,numChannels);
        
        for cc = 1:numChannels
        
            [m,n,o]=size(rawStack{cc}); %M is the original matrix
            
            m = binning.*floor(m/binning);
            n = binning.*floor(n/binning);
            o = binning.*floor(o/binning);
            
            % dimension 1 binning
            
            xyBinStack = zeros(m/binning,n/binning,o/binning);
            
            for ll = 1:o
                
                ll
            
                container = squeeze(rawStack{cc}(1:m,1:n,ll));
                
                container = sum(reshape(container,binning,[]),1);
                container = reshape(container,m/binning,[]);
                container = permute(container,[2,1]);
                
                container = sum(reshape(container,binning,[]),1);
                container = reshape(container,n/binning,[]);

                xyBinStack(:,:,ll) = permute(container,[2,1]);
                
            end
            
            xyBinStack = permute(xyBinStack,[3,2,1]);
            xyBinStack = sum(reshape(xyBinStack,binning,[],m/binning),1);
            xyBinStack = reshape(xyBinStack,o/binning,[],m/binning);
            xyBinStack = permute(xyBinStack,[3,2,1]);
            
            binnedStack{cc} = xyBinStack./binning.^3;
            
        end

    else
        
        binnedStack = rawStack;
        
    end
    
    clear('rawStack')
    
    [stackSizeY,stackSizeX,stackSizeZ] = size(binnedStack{segChannel});
    voxelSizeY = rawVoxelSizeY.*binning;
    voxelSizeX = rawVoxelSizeY.*binning;
    voxelSizeZ = rawVoxelSizeZ.*binning;
    voxelVol = voxelSizeY.*voxelSizeX.*voxelSizeZ;
    
    maxProjCellYY{kk} = squeeze(max(binnedStack{segChannel},[],1));
    maxProjCellXX{kk} = squeeze(max(binnedStack{segChannel},[],2));
    maxProjCellZZ{kk} = squeeze(max(binnedStack{segChannel},[],3));
    
    volCell = cell(1,threshSteps);
    centrCell = cell(1,threshSteps);
    vxlIdxCell = cell(1,threshSteps);
    coordsVxlIdxCell = cell(1,threshSteps);
    imagesCell = cell(1,threshSteps);
    
    numRegVec = zeros(1,threshSteps);
    
    featureVecCell = cell(1,threshSteps);
    
    for tt = 1:threshSteps
        
        disp(threshSteps - tt);
        
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
        
    end
    
    % --- Linking of objects
    
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
    
    
    % -- Quantification based on extracted regions
    
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
        nucCent{kk}(ll) = ...
            nucCent{kk}(l).*[voxelSizeY,voxelSizeX,voxelSizeZ];
    end
    nucBBox{kk} = nucleiBoundingBox;
    nucVxlIdx{kk} = nucleiVxlIdx;
    nucImg{kk} = nucleiImage;
    
    nucInt{kk} = cell(1,numChannels);
    cytoInt{kk} = cell(1,numChannels);
    
    for cc = 1:numChannels
        
        % Extract further properties only for the valid regions
        nucleiProps = regionprops(nucleiRegions,binnedStack{cc},...
            'MeanIntensity');
        
        nucleiMeanInt = [nucleiProps(:).MeanIntensity];
        nucInt{kk}{cc} = nucleiMeanInt;
        
        cytoInt{kk}{cc} = zeros(1,numNuc(kk));
        
    end
    
    
    % Make masks to determine cytoplasmic intensity
    
    totalDil = dilateSpacer+dilateMeasure;
    
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
        extMaxY = min(rawStackSizeY,fullExtMaxY);
        yHiDiff = fullExtMaxY - extMaxY;
        extMaxX = min(rawStackSizeX,fullExtMaxX);
        xHiDiff = fullExtMaxX - extMaxX;
        extMaxZ = min(rawStackSizeZ,fullExtMaxZ);
        zHiDiff = fullExtMaxZ - extMaxZ;
        
        % Extended bounding box size
        extSizeY = extMaxY - extMinY + 1;
        extSizeX = extMaxX - extMinX + 1;
        extSizeZ = extMaxZ - extMinZ + 1;
        
        % Inclusion mask
        inclMask = zeros(extSizeY,extSizeX,extSizeZ);
        inclMask((1+totalDil-yLoDiff):(end-totalDil+yHiDiff),...
            (1+totalDil-xLoDiff):(end-totalDil+xHiDiff),...
            (1+totalDil-zLoDiff):(end-totalDil+zHiDiff))...
            = nucImg{kk}{nn};
        
        % Exclusion mask
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
            
            cytoInt{kk}{cc}(nn) = mean(cutoutImage(measureMask));
            
        end
        
    end
    
    
    % --- Nearest neighbor distances
                
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
    
    if ~plot_nearest_neighbor
        
        % --- plotting of nearest neighbor distances
        
        figure(2)
                
        plot(NN_median(kk).*[1,1],[0,1],'r--')
        
        hold on
        
        plot(sort(distances),linspace(0,1,numel(distances)),'k-o')
        
        hold off
        
        xlabel('Individual nucleus NN distance [\mum]')
        ylabel('Cumulative probability')
        
    end
    
    
    % Upon successful completion, set error flag to false, no error!
    errorFlagVec(kk) = false;
    
    %     catch message
    %
    %         errorFlagVec(kk) = true;
    %         errorMessages{kk} = message;
    %
    %         fprintf('Image segmentation error in step %d.\n',kk)
    %
    %     end
    
end


%% --- plot max projection

figure(1)

subplot(2,3,1)

imagesc(maxProjCellYY{1}.')

colormap(gray)

hold on

set(gca,'XLim',[0,stackSizeX],'YLim',[0,stackSizeZ],'YDir','normal')

xlabel('x [\mum]')
ylabel('z [\mum]')

hold off

daspect([voxelSizeZ,voxelSizeY,1])

title('Y maximum projection')


subplot(2,3,2)

imagesc(maxProjCellXX{1}.')

colormap(gray)

hold on

set(gca,'XLim',[0,stackSizeY],'YLim',[0,stackSizeZ],'YDir','normal')

xlabel('y [\mum]')
ylabel('z [\mum]')

hold off

title('X maximum projection')

daspect([voxelSizeZ,voxelSizeY,1])



subplot(2,3,3)

imagesc(maxProjCellZZ{1})

colormap(gray)

hold on

set(gca,'XLim',[0,stackSizeX],'YLim',[0,stackSizeY])

xlabel('x [\mum]')
ylabel('y [\mum]')

hold off

title('Z maximum projection')

daspect([voxelSizeY,voxelSizeX,1])








subplot(2,3,4)

imagesc(maxProjCellYY{1}.')

colormap(gray)

hold on

set(gca,'XLim',[0,stackSizeX],'YLim',[0,stackSizeZ],'YDir','normal')

xlabel('x [\mum]')
ylabel('z [\mum]')

plot(yCoords./voxelSizeY,zCoords./voxelSizeZ,'r+')

hold off

daspect([voxelSizeZ,voxelSizeY,1])

title('Y maximum projection')


subplot(2,3,5)

imagesc(maxProjCellXX{1}.')

colormap(gray)

hold on

set(gca,'XLim',[0,stackSizeY],'YLim',[0,stackSizeZ],'YDir','normal')

xlabel('y [\mum]')
ylabel('z [\mum]')

plot(xCoords./voxelSizeX,zCoords./voxelSizeZ,'r+')

hold off

title('X maximum projection')

daspect([voxelSizeZ,voxelSizeY,1])



subplot(2,3,6)

imagesc(maxProjCellZZ{1})

colormap(gray)

hold on

set(gca,'XLim',[0,stackSizeX],'YLim',[0,stackSizeY])

xlabel('x [\mum]')
ylabel('y [\mum]')

plot(yCoords./voxelSizeY,xCoords./voxelSizeX,'r+')

hold off

title('Z maximum projection')

daspect([voxelSizeY,voxelSizeX,1])
