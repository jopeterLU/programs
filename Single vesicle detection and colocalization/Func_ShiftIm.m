function [I_DC,Mask,MCC,PCC] = Func_ShiftIm(im,Mask,side,offset)

cInd_subIm = size(im,2)/2;
rInd_subIm = size(im,1);

redIm = im(:,1:cInd_subIm);
blueIm = im(:,cInd_subIm+1:end);

redMask = Mask(:,1:cInd_subIm);
blueMask = Mask(:,cInd_subIm+1:end);

if  strcmp(side,'left')==1
    HalfMask = redMask;
    UHalfMask = blueMask;
    rIm = redIm;
    fIm = blueIm;
elseif strcmp(side,'right')==1
    HalfMask = blueMask;
    UHalfMask = redMask;
    rIm = blueIm;
    fIm = redIm;
end


% Shift matrix
Im_DriftCorr = zeros(size(HalfMask));
if offset(1)<0
    st = (-1)*offset(1);
    Im_DriftCorr(1:end-st,:) = fIm(st+1:end,:);
    UHalfMask(1:end-st,:) = UHalfMask(st+1:end,:);
else
    st = offset(1);
    Im_DriftCorr(st+1:end,:) = fIm(1:end-st,:);
    UHalfMask(st+1:end,:) = UHalfMask(1:end-st,:);
end
if offset(2)<0
    st = (-1)*offset(2);
    Im_DriftCorr(:,1:end-st) = Im_DriftCorr(:,st+1:end);
    UHalfMask(:,1:end-st) = UHalfMask(:,st+1:end);
else
    st = offset(2);
    Im_DriftCorr(:,st+1:end) = Im_DriftCorr(:,1:end-st);
    UHalfMask(:,st+1:end) = UHalfMask(:,1:end-st);
end

if  strcmp(side,'left')==1
    I_DC = [rIm Im_DriftCorr];
    Mask = [HalfMask UHalfMask];
    
elseif strcmp(side,'right')==1
    I_DC = [Im_DriftCorr rIm];
    Mask = [UHalfMask HalfMask];
end

MCC = sum(sum(rIm.*Im_DriftCorr))/(sqrt(sum(sum(rIm.^2))*sum(sum(Im_DriftCorr.^2))));
%PCC = sum(sum((rIm-mean2(rIm)).*(Im_DriftCorr-mean2(Im_DriftCorr))))/sqrt(sum(sum((rIm-mean2(rIm)).^2))*sum(sum((Im_DriftCorr-mean2(Im_DriftCorr)).^2)));
PCC = corr2(rIm,Im_DriftCorr);


end

