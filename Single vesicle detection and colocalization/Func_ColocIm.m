function [ CoStruct,N_co ] = Func_ColocIm(r_shift,c_shift,redSpots,blueSpots,RedData,BlueData,minS)

N_co = 0;
CoStruct = [];

for l = 1:redSpots
        r_cent = RedData(l).rCent;
        c_cent = RedData(l).cCent; 

        lim = [r_cent-r_shift,r_cent+r_shift,c_cent-c_shift,c_cent+c_shift];
        nv = 0;
        for o = 1:blueSpots
            r = BlueData(o).rCent;
            c = BlueData(o).cCent;   

             if  lim(1) <= r && r <= lim(2) && lim(3) <= c && c <= lim(4) && nv ==0....
                 RedData(l).size > minS && BlueData(o).size > minS
                N_co = N_co + 1;
                adding.redTot = RedData(l).I_tot;
                adding.redMax = RedData(l).Imax;
                adding.redSize = RedData(l).size;
                adding.redRow = r_cent;
                adding.redCol = c_cent;
                adding.blueTot = BlueData(o).I_tot;
                adding.blueMax = BlueData(o).Imax;
                adding.blueSize = BlueData(o).size;
                adding.blueRow = r;
                adding.blueCol = c;
                adding.redim = RedData(l).slide;
                adding.blueim = BlueData(o).slide;
                CoStruct = [CoStruct adding];  
                clear adding;
                clear BlueData(o);
                nv = 1;
             end      
        end
        
end

end

