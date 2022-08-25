%% Script for a folder with separate images

% Input: .tif stack file 
% Output: .mat-file containing structures with information about all
% detected spots

% Copyright: Alexandra Andersson
% Lund University, Division of Physical Chemistry
% 220823


%% Define parameters

radius = 100; % radius of rolling ball
dx = 2; % spacing RB
SNR_red = 3;
SNR_blue = 3;
minSize = 25;
maxSize = 500;
rel_max = 0.1;
pixA = 0.0469;
r_shift = 5;
c_shift = 5;
th_bleedthrough = 0;

imRow = 1024;
imCol = 1024;
nChannel = 2;
chanRow = 1024;
chanCol = 512;
cropUp = 200; cropLow = 200; cropLeft = 100; cropRight = 100;


prompt = {'Full image rows','Full image columns','Number of channels',...
    'Channel size rows','Channel size columns','cropUp','cropLow','cropLeft',...
    'cropRight','pixA','x_shift','y_shift','threshold bleedthrogh'};
name ='Choose image type';
numlines = 13;
options.Resize = 'on';
size_wind = [1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50; 1 50;];
defaultanswer = {mat2str(imRow),mat2str(imCol),mat2str(nChannel),mat2str(chanRow),...
    mat2str(chanCol),mat2str(cropUp),mat2str(cropLow),mat2str(cropLeft),...
    mat2str(cropRight),mat2str(pixA),mat2str(c_shift),mat2str(r_shift),mat2str(th_bleedthrough)};
upd = inputdlg(prompt,name,size_wind,defaultanswer,options);

imRow = str2double(upd{1});
imCol= str2double(upd{2});
nChannel = str2num(upd{3});
chanRow = str2double(upd{4});
chanCol = str2num(upd{5});
cropUp = str2num(upd{6});
cropLow = str2num(upd{7});
cropLeft = str2num(upd{8});
cropRight = str2num(upd{9});
pixA = str2num(upd{10});
c_shift = str2num(upd{11});
r_shift = str2num(upd{12});
th_bleedthrough = str2num(upd{13});

close

rnd = linspace(cropUp,chanRow-cropLow,chanRow-cropUp-cropLow+1);
cnd = linspace(cropLeft,chanCol-cropRight,chanCol-cropLeft-cropRight+1);
f_cutoff = 500;

list = {'Yes,left side','Yes, right side','No','Yes, both sides'};
[ICidx,tf] = listdlg('PromptString',{'Intensity copmpared colocalization?',''},...
    'SelectionMode','single','ListString',list);

DC = questdlg('Drift-correct?',...
    'Yes','No');
if strcmp(DC,'No')==1
    prompt = {'x','y'};
    name ='Offset';
    numlines = 2;
    options.Resize = 'on';
    size_wind = [1 50;1 50;];
    defaultanswer = {mat2str(0),mat2str(0)};
    DCupd = inputdlg(prompt,name,size_wind,defaultanswer,options);
    offset = [str2double(DCupd{1}) str2double(DCupd{2})];
end

%% Load stack
 
[FileName,FilePath]=uigetfile('.tif');

% For stack image
ExPath = [FilePath FileName];
fileinfo = imfinfo(ExPath);
f_nr = numel(fileinfo);

%% Detect spots and colocalization


CoData = cell(f_nr,1); % Struct containing first info of all spots codetected for whole stack
RedData = cell(f_nr,1); % Struct containing first info of all red spots detected for whole stack
BlueData = cell(f_nr,1); % Struct containing first info of all blue spots detected for whole stack

red_bg_spot = cell(f_nr,1);
blue_bg_spot = cell(f_nr,1);

redSpots = zeros(1,f_nr);
blueSpots = zeros(1,f_nr);

Im_Coloc = zeros(1,f_nr);
Im_colocInd = zeros(1,f_nr);

MCC = zeros(1,f_nr);
MCC_NoDC = zeros(1,f_nr);
PCC = zeros(1,f_nr);
PCC_ML = zeros(1,f_nr);
%%   
cont = 'Yes';


for k = 1:f_nr

    I_org = double(imread(ExPath, k, 'Info', fileinfo));

    I_LP = Func_ILPF(I_org,f_cutoff);
    BG = Func_RollingBall(I_LP,radius,dx);
    I_BGsub = I_LP-BG;
    Im  = Func_BGcorr(I_BGsub);

    % Asked about SNR update for every image 
    if strcmp(cont,'Yes')==1
    [Mask,rel_max,SNR_red,SNR_blue] = Func_MaskInteractive(Im,SNR_blue,SNR_red,rnd,cnd,minSize,maxSize,chanRow,chanCol,nChannel,rel_max);
    % Use non-interactive Mask 
    elseif strcmp(cont,'No')==1
    [Mask] = Func_MaskNonInt(Im,SNR_blue,SNR_red,rnd,cnd,minSize,maxSize,chanRow,chanCol,nChannel);
    end 
%         im = Im; mask = Mask;
    if strcmp(DC,'Yes')==1
        [Im,Mask,offset,MCC(k),PCC(k)] = Func_DriftCorr(Im,Mask,'left');
    elseif strcmp(DC,'No')==1
        [Im,Mask,MCC(k),PCC(k)] = Func_ShiftIm(Im,Mask,'left',offset);
    end

    MCC_NoDC(k) = sum(sum(Im(:,1:chanCol).*Im(:,chanCol+1:end)))/(sqrt(sum(sum(Im(:,1:chanCol).^2))*sum(sum(Im(:,chanCol+1:end).^2))));
    [red_bg_spot{k},blue_bg_spot{k},c] = Func_BG_spot(Mask,Im,chanCol,chanRow,minSize,maxSize);

    % Detect ves one channel,calculate intensity in second channel    
    if not (ICidx ==3)

        if isequal (ICidx,1)
             [leftSpots(k),leftData{k},leftCo{k},leftNonCo{k},leftCoInd(k)] = Func_IntCalcVes(Im(:,1:chanCol),Mask(:,1:chanCol),rnd,cnd,k,Im(:,chanCol+1:end),Mask(:,chanCol+1:end),red_bg_spot{k},blue_bg_spot{k},c,maxSize,minSize,SNR_blue);                
        elseif isequal(ICidx,2)
            [rightSpots(k),rightData{k},rightCo{k},rightNonCo{k},rightCoInd(k)] = Func_IntCalcVes(Im(:,chanCol+1:end),Mask(:,chanCol+1:end),rnd,cnd,k,Im(:,1:chanCol),Mask(:,1:chanCol),blue_bg_spot{k},red_bg_spot{k},c,maxSize,minSize,SNR_red); 
        elseif isequal(ICidx,4)
            [leftSpots(k),leftData{k},leftCo{k},leftNonCo{k},leftCoInd(k)] = Func_IntCalcVes(Im(:,1:chanCol),Mask(:,1:chanCol),rnd,cnd,k,Im(:,chanCol+1:end),Mask(:,chanCol+1:end),red_bg_spot{k},blue_bg_spot{k},c,maxSize,minSize,SNR_blue);                
            [rightSpots(k),rightData{k},rightCo{k},rightNonCo{k},rightCoInd(k)] = Func_IntCalcVes(Im(:,chanCol+1:end),Mask(:,chanCol+1:end),rnd,cnd,k,Im(:,1:chanCol),Mask(:,1:chanCol),blue_bg_spot{k},red_bg_spot{k},c,maxSize,minSize,2*SNR_red);               
        end
    end

    % Detect vesicles in separate channels
    [redSpots(k),redBG(k),redSTD(k),RedStruct] = Func_IntCalc(Im(:,1:chanCol),Mask(:,1:chanCol),rnd,cnd,k,pixA,th_bleedthrough);
    [blueSpots(k),blueBG(k),blueSTD(k),BlueStruct] = Func_IntCalc(Im(:,chanCol+1:end),Mask(:,chanCol+1:end),rnd,cnd,k,pixA,0);

    RedData{k} = RedStruct;
    BlueData{k} = BlueStruct;

    % Channel-dependent colocalization within image

    [CoStruct,Im_Coloc(k)] = Func_ColocIm(r_shift,c_shift,redSpots(k),blueSpots(k),RedStruct,BlueStruct,0*pixA);
    CoData{k} = CoStruct;
    Im_colocInd(k) = 200.*Im_Coloc(k)./(redSpots(k)+blueSpots(k));
    RB = redSpots./blueSpots;

    if isempty(CoStruct)==1 

        if strcmp(cont,'Yes')==1
            answer = questdlg('No colocalized vesicles found. Continue?');
            if strcmp(answer,'No')==1
                break
            end
            close
            cont = questdlg('Continue interactively?');
            close 
        end

    elseif isempty(CoStruct)==0

        redIm = Mask(rnd,cnd);
        blueIm = Mask(rnd,chanCol+1+cnd);
        Im_c = imfuse(redIm,blueIm,'falsecolor','ColorChannels',[1 0 2]);

        if strcmp(cont,'Yes')==1
            % Plotting coloc of specific image
            admin = min(min(Im));
            admax = max(max(Im));
            admean = mean2(Im);
            clims = [admin admean];
            figure('units','normalized','outerposition',[0 0 1 1])
            imagesc(Im_c,clims)
            hold on

            for w = 1:length(CoStruct)
                viscircles([CoStruct(w).redCol-cropLeft,CoStruct(w).redRow-cropUp],10,'Color','w')
            end 
            answer = questdlg('Colocalized vesicles found. Continue?');
            if strcmp(answer,'No')==1
                break
            end
            close
            cont = questdlg('Continue interactively?');
            close       
        end
    end

    disp([mat2str(k),' of ', mat2str(f_nr)])
end


%% Save matfile

[file,path] = uiputfile('filename.mat');       

disp('Saving...')
save('-v7.3',[path '\' file]);
disp('Saved.')


