close all
% Tiff-stack to analyze (cropped with background subtracted)
% Displays a dialog box from where to choose the images to be analyzed
[FileName,PathName] = uigetfile('*.tif','Select the image-file');
if FileName==0
    return
end
filename=strcat(PathName,FileName);

% Parameters affecting the detection of spots, please adjust to detect all
% spots
cell_max=6; % diameter of a spot in pixels
I_thr=50; % threshold intensity

% The spots are detected from the first image in the stack
I=double(imread(filename,1));

I=double(I);
h=fspecial('gaussian',5,2); % Gaussian filter to reduce noise
I2=imfilter(I,h,'replicate');

out=pkfnd(I2,I_thr,cell_max);   % Detection of spots using pkfnd

% This figure shows the detected spots, please adjust cell_max and I_thr
% until the detection is satisfactory
figure(1)
imshow(I2,[0 I_thr])
viscircles(out,ones(1,size(out,1))*cell_max,'LineWidth',0.25,'Color','r');

num_spots=length(out);
x=[1:size(I,2)];
y=[1:size(I,1)];
[Xx,Yy]=meshgrid(x,y);

R2=cell_max^2;
BW=zeros(size(I));
for i=1:num_spots
   BW(((Xx-out(i,1)).^2+(Yy-out(i,2)).^2)<R2)=1;
end

% BW is an image with a circle (radius cell_max) around each detected spot
figure(2)
imshow(BW,[])

BW2=imdilate(BW,strel('disk',cell_max));
BW3=imdilate(BW,strel('disk',2*cell_max));
BW2=BW3.*(1-BW2);   

% BW2 is a ring around each detected spot to obtain the local backround from
figure(3)
imshow(BW2,[])

% Ispot is the total intensity from a spot
Ispot=(sum(I(BW==1))-mean(I(BW2==1))*sum(sum(BW==1)))/num_spots
