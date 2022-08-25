function [Mask,rel_max,SNR_red,SNR_blue] = Func_MaskInteractive(im,SNR_blue,SNR_red,rnd,cnd,minSize,maxSize,chanRow,chanCol,nChannel,rel_max)

% Function called to update threshold of filtered image

% In: Image to be analyzed, SNR, indices, max and min sizes, channel sizes,
% and contrast adjustment factor

% Out: Image with only spots (Mask), updated contrast adjustment and SNRs

if nChannel == 2
    
    cInd_subIm = chanCol;
    rInd_subIm = chanRow;

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

    answer = 'Yes';
    admin = min(min(im));
    admax = max(max(im));
    admean = mean2(im);
    while strcmp(answer,'Yes')==1
        
        blueTh = SNR_blue.*2.*repmat(BlueNoise,[rInd_subIm,1]);
        redTh = SNR_red.*2.*repmat(RedNoise,[rInd_subIm,1]);

        th_mask = [redTh blueTh];
        Mask = (im > th_mask);
        Mask(:,1:cInd_subIm) = Mask(:,1:cInd_subIm).*roiMask;
        Mask(:,cInd_subIm+1:end) = Mask(:,cInd_subIm+1:end).*roiMask;
        enlarge = ones(round(5));
        Mask = imdilate(Mask,enlarge);
        
        SpotInf = bwconncomp(Mask);
 
        for i = 1:SpotInf.NumObjects
            indList = length(SpotInf.PixelIdxList{i});

            if indList < minSize | indList > maxSize
                Mask(SpotInf.PixelIdxList{i}) = 0;
            end
        end

        clims = [admin rel_max.*admax];
        [B,~] = bwboundaries(Mask);
        figure('units','normalized','outerposition',[0 0 1 1])
        imagesc(im,clims)
        colormap(jet)
        cb = colorbar();
        cb.Ruler.Scale = 'log';
        cb.Ruler.MinorTick = 'on';

        hold on
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5)
        end

        prompt = {'SNR red','SNR blue','Relative max scale'};
        name ='Change SNR';
        numlines = 3;
        defaultanswer = {mat2str(SNR_red),mat2str(SNR_blue),mat2str(rel_max)};
        upd = newid(prompt,name,numlines,defaultanswer);
        answer = questdlg('Update?','Update SNR','No');
        if strcmp(answer,'Yes')==1
            SNR_red = str2double(upd{1});
            SNR_blue = str2double(upd{2});
            rel_max = str2double(upd{3});

        end
        close
        if isempty(sum(Mask))
            out = [];
            disp('No spots found')
            return;
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

    answer = 'Yes';
    admin = min(min(im));
    admax = max(max(im));
    admean = mean2(im);

    while strcmp(answer,'Yes')==1

        Th = SNR_red.*2.*repmat(Noise,[imRow,1]);
        Mask = (im > Th);
        Mask = Mask.*roiMask;
        enlarge = ones(round(5));
        Mask = imdilate(Mask,enlarge);

        clims = [admin rel_max.*admax];
        [B,~] = bwboundaries(Mask);
        figure('units','normalized','outerposition',[0 0 1 1])
        imagesc(im,clims)
        colormap(jet)
        cb = colorbar();
        cb.Ruler.Scale = 'log';
        cb.Ruler.MinorTick = 'on';

        hold on
        for k = 1:length(B)
           boundary = B{k};
           plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 1.5)
        end

        prompt = {'SNR','Relative max scale'};
        name ='Change SNR';
        numlines = 2;
        defaultanswer = {mat2str(SNR_red),mat2str(rel_max)};
        upd = newid(prompt,name,numlines,defaultanswer);
        answer = questdlg('Update?','Update SNR','No');
        if strcmp(answer,'Yes')==1
            SNR_red = str2double(upd{1});
            rel_max = str2double(upd{2});

        end
        close
        if isempty(sum(Mask))
            out = [];
            disp('No spots found')
            return;
        end
    end

    SpotInf = bwconncomp(Mask);
    for i = 1:SpotInf.NumObjects
        indList = length(SpotInf.PixelIdxList{i});

        if indList < minSize | indList > maxSize
            Mask(SpotInf.PixelIdxList{i}) = 0;
        end
    end

end
close all
out = Mask;
