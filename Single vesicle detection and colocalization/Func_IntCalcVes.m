function [spots,SpotStatStruct,CoStruct,NonCoStruct,n_co] = Func_IntCalcVes(im,Mask,rnd,cnd,slide,im2,Mask2,bg_ch1,bg_ch2,c,maxSize,minSize,SNR)


%%
bg_ch1 = mean(bg_ch1,2);
bg_ch2 = mean(bg_ch2,2);
bgPch1 = polyfit(c(:,1),bg_ch1,1);
bgPch2 = polyfit(c(:,1),bg_ch2,1);

x_size = linspace(1,maxSize,maxSize);
int_ch1 = bgPch1(1)*x_size+bgPch1(2);
int_ch2 = bgPch2(1)*x_size+bgPch2(2);

roiMask = zeros(size(im,1),size(im,2));
roiMask(rnd,cnd) = 1;
SpotStatStruct = [];
CoStruct = [];
NonCoStruct = [];
co_count = 0;

%% Masking

Mask = Mask.*roiMask;
Mask2 = Mask2.*roiMask;
InvMask = ~Mask;
InvIm = im.*InvMask;
InvMask2 = ~Mask2;
InvIm2 = im2.*InvMask2;

BG_mean_ch1 = mean2(InvIm);
BG_mean_ch2 = mean2(InvIm2);
BG_std_ch1 = std2(InvIm)/(sqrt(size(im,1)*size(im,2)));

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
        I_p = zeros(1,pixL);
        for j = 1:pixL
            I(j) = im(rInd(j),cInd(j));
            I_p(j) = im2(rInd(j),cInd(j));
        end
        I_tot = sum(I);
        [I_max,~] = max(I);  
        
        Ip_tot = sum(I_p);
        [Ip_max,~] = max(I_p);  
        
        adding.Imax_ch1 = I_max-BG_mean_ch1;
        adding.Imax_ch2 = Ip_max-BG_mean_ch2;
        adding.ch2_ch1_max = (Ip_max-BG_mean_ch2)/(I_max-BG_mean_ch1);
        adding.I_tot_ch1 = I_tot-int_ch1(pixL);
        adding.I_tot_ch2 = Ip_tot-int_ch2(pixL);       
        adding.ch2_ch1_tot = (Ip_tot-int_ch2(pixL))/(I_tot-int_ch1(pixL));
        adding.slide = slide;
        adding.rInd = rInd;
        adding.cInd = cInd;
        adding.rCent = rCent;
        adding.cCent = cCent;
        adding.size = pixL*0.0469;
        
        if Ip_tot-int_ch2(pixL) > int_ch2(minSize)*SNR && (Ip_tot-int_ch2(pixL))/(I_tot-int_ch1(pixL))>0.1
            co_count = co_count+1;
            adding.co = 'Yes';
            CoStruct = [CoStruct adding];
        else
            adding.co = 'No';
            NonCoStruct = [NonCoStruct adding];
        end

        SpotStatStruct = [SpotStatStruct adding];
        clear adding

    end
    spots = length(SpotStatStruct);
    n_co = 100*co_count/spots;

end
