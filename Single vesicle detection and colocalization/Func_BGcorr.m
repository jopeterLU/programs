function [ BGsub ] = Func_BGcorr( BGsub )

for i = 1:size(BGsub,1)
    for j = 1:size(BGsub,2)
        if BGsub(i,j) < 0
            BGsub(i,j) = 1;
           
        elseif BGsub(i,j) == 0
            BGsub(i,j) = 1;
        end 
    end
end      

BGsub = imgaussfilt(BGsub);
end

