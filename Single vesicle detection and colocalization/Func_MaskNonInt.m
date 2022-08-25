function Mask = Func_MaskNonInt(im,SNR_blue,SNR_red,rnd,cnd,minSize,maxSize,imRow,imCol,nChannel)

% Function called in GUI to update threshold of filtered image

% In: Image to be analyzed, factor to multiply original threshold with,
% indices, min and max size, cropped area

% Out: Image with only spots (Mask)

if nChannel == 2
    
    cInd_subIm = imCol;
    rInd_subIm = imRow;

    redIm = im(:,1:cInd_subIm);
    blueIm = im(:,cInd_subIm+1:end);

    RedMed = median(redIm);
    BlueMed = median(blueIm);

    RedMin = min(redIm);
    BlueMin = min(blueIm);

    RedMax = max(redIm);
    BlueMax = max(blueIm);

    RedNoise = RedMed-RedMin;
    BlueNoise = BlueMed-BlueMin;

    % Choose ROI in Mask
    roiMask = zeros(rInd_subIm,cInd_subIm);
    roiMask(rnd,cnd) = 1;

    blueTh = SNR_blue.*2.*repmat(BlueNoise,[rInd_subIm,1]);
    redTh = SNR_red.*2.*repmat(RedNoise,[rInd_subIm,1]);

    th_mask = [redTh blueTh];
    Mask = (im > th_mask);
    Mask(:,1:cInd_subIm) = Mask(:,1:cInd_subIm).*roiMask;
    Mask(:,cInd_subIm+1:end) = Mask(:,cInd_subIm+1:end).*roiMask;
    enlarge = ones(round(5));
    Mask = imdilate(Mask,enlarge);

    if isempty(sum(Mask))
        out = [];
        disp('No spots found')
        return;
    end

    SpotInf = bwconncomp(Mask);
    for i = 1:SpotInf.NumObjects
        indList = length(SpotInf.PixelIdxList{i});

        if indList < minSize | indList > maxSize
            Mask(SpotInf.PixelIdxList{i}) = 0;
        end
    end

elseif nChannel ==1
    
    Med = median(im);
    Min = min(im);
    Max = max(im);
    Noise = Med-Min;

    % Choose ROI in Mask
    roiMask = zeros(imRow,imCol);
    roiMask(rnd,cnd) = 1;

    Th = SNR_red.*2.*repmat(Noise,[imRow,1]);
    Mask = (im > Th);
    Mask = Mask.*roiMask;
    enlarge = ones(round(5));
    Mask = imdilate(Mask,enlarge);

    if isempty(sum(Mask))
        out = [];
        disp('No spots found')
        return;
    end

    SpotInf = bwconncomp(Mask);
    for i = 1:SpotInf.NumObjects
        indList = length(SpotInf.PixelIdxList{i});

        if indList < minSize | indList > maxSize
            Mask(SpotInf.PixelIdxList{i}) = 0;
        end
    end

end
end
