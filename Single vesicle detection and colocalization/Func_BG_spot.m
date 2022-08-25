function [red_bg_spot,blue_bg_spot,c] = Func_BG_spot(Mask,Im,chanCol,chanRow,minSize,maxSize)
%Check spot-BG in Mask
        redMask = not(Mask(:,1:chanCol));
        blueMask = not(Mask(:,chanCol+1:end));
        redIm = Im(:,1:chanCol);
        blueIm = Im(:,chanCol+1:end);
        redCI = redMask.*redIm; blueCI = blueMask.*blueIm;
        spotSize = round(sqrt(linspace(minSize,maxSize,10)));
        red_bg_spot = zeros(length(spotSize),10);
        blue_bg_spot = zeros(length(spotSize),10);
        c = zeros(length(spotSize),10);
        for p = 1:length(spotSize)
            for q = 1:10
                c_idx = randsample(1:chanCol,spotSize(p),false);
                r_idx = randsample(1:chanRow,spotSize(p),false);
                for r = 1:length(c_idx)
                    for s = 1:length(r_idx)
                        red_bg_spot(p,q) = red_bg_spot(p,q)+(redCI(r_idx(s),c_idx(r)));
                        blue_bg_spot(p,q) = blue_bg_spot(p,q)+(blueCI(r_idx(s),c_idx(r)));
                        c(p,q) = c(p,q)+1;
                    end
                end
            end
        end
            
end

