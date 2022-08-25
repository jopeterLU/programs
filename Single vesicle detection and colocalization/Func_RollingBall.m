function RB_smooth = Func_RollingBall(Im,radius,dx)
% Rolling ball filter for ball below input
% Im is image where rolling ball algorithm should be performed
% radius is the radius of the ball rolling under graph
% dx is spacing of sample points

N_col = size(Im,2);
N_row = size(Im,1);
RB_out = zeros(N_row,N_col); % Define BG that will be output
K = floor(radius/dx);   %number of neghibours

for i = 1:N_col
    Im_c = Im(:,i);   
    BG_col = Im_c;  % Start at the value  
    
    for k = 1:K %Constrain from kth left/right neighbors
      st = k + 1; %start marker
      en = N_row - k;  %end marker
      V = Im_c-sqrt(radius^2 - (k*dx)^2);
      BG_col(1:en) = min(  BG_col(1:en),V(st:N_row) );  %left
      BG_col(st:N_col) = min( BG_col(st:N_row),V(1:en) );  %right
    end
    
    BG_line = BG_col + radius; % Add so that bg level lifts from middle of ball to below graph
    RB_out(:,i) = BG_line;
%     RB_smooth(:,i) = smooth(1:length(BG_line),BG_line,0.5,'rloess'); 
end
% 
% K = 0.045*ones(5);
% RB_smooth = conv2(RB_out,K,'same');
 RB_smooth = medfilt2(RB_out);
end

