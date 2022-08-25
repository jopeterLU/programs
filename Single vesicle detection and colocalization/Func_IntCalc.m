function [spots,BG_mean,BG_std,SpotStatStruct] = Func_IntCalc(im,Mask,rnd,cnd,slide,pixA,th_bleedthrough)

% Function that calculates calculates mean BG intensity and intensities in
% the spots
% In: im, Mask, indices, slide, pixel area, threshold bleedthrough
% Out: number of spots; BGmean, BGstd; structure over spots (I, ind, cent, size, slide)

%%

roiMask = zeros(size(im,1),size(im,2));
roiMask(rnd,cnd) = 1;
SpotStatStruct = [];

%% Masking

Mask = Mask.*roiMask;
InvMask = ~Mask;
InvIm = im.*InvMask;

BG_mean = mean2(InvIm);
BG_std = std2(InvIm)/(sqrt(size(im,1)*size(im,2)));

MaskInf = bwconncomp(Mask);

% Find median value for each spot
    for i = 1:MaskInf.NumObjects

        pixL = length(MaskInf.PixelIdxList{i});
        Pix = MaskInf.PixelIdxList{i};
        rInd = mod(Pix,size(Mask,1));
        cInd = floor(Pix/size(Mask,1))+1;
        rCent = median(rInd);
        cCent = median(cInd);

        I = zeros(1,pixL);

        for j = 1:pixL
            I(j) = im(rInd(j),cInd(j));
        end
        I_tot = sum(I);
        [I_max,~] = max(I);  
        
        if I_max > th_bleedthrough
            adding.Imax = I_max;
            adding.I_tot = I_tot;
            adding.slide = slide;
            adding.rInd = rInd;
            adding.cInd = cInd;
            adding.rCent = rCent;
            adding.cCent = cCent;
            adding.size = pixL*pixA;
            SpotStatStruct = [SpotStatStruct adding];
            clear adding
        end

    end
    spots = length(SpotStatStruct);

end
