%%This script is how we define the training data and develop the ROC curve
close all
clear all
clc

addpath('~/Documents/GitHub/UniversalCode/');
addpath('~/Documents/GitHub/Islet_Heterogeneity/')
%chose islets to analyze
filename = ["three","sample","five", "two","one"];
csvname = ["H2BmCherry Ucn3GCaMP-3_Detailed.csv", "H2BmCherry Ucn3GCaMP sample_Detailed.csv", "H2BmCherry Ucn3GCaMP-5_Detailed.csv", "H2BmCherry Ucn3GCaMP-2_Detailed.csv", "H2BmCherry Ucn3GCaMP-1_Detailed.csv"]
datapath = ['/Volumes/Briggs_10TB/Merrin/Confocal/'] 
savepath = ['~/Documents/GitHub/ST_Analysis/Data/']


mm = 1; %set which type of mask
masktypes = {'Bad', 'Medium', 'Good'}
masktype = masktypes{mm}
%initialize which time index to use nuclear locations
%timetouse = 243%three: 217 - somewhat wiggly. %sample: 194; %five: 243%two: 363; %one: 193
timetouse = [217, 194, 243, 363, 193];
%Number of cells per islet for the training case
perislet = 6
%set seed

for kt = 1:length(filename)
ca_files = dir(strrep(strjoin([datapath filename(kt) '/' '*C1*.tif']),' ',''));
nuc_files = dir(strrep(strjoin([datapath filename(kt) '/' '*C2*.tif']),' ',''));

    %% Import Images %%
    %import calcium
    for i = 1:length(ca_files)
        fulldatapath = strrep(strjoin([datapath filename(kt) '/' ca_files(i).name]),' ','');
        image = imread(fulldatapath);
        Islet_vid(:,:,i) = image; %calcium files are loaded as a X pixel x Y pixel x Time
        %imshow(imcomplement(image));, drawnow
        fullnucpath = strrep(strjoin([datapath filename(kt) '/' nuc_files(i).name]),' ', '');
        nuimage = imread(fullnucpath);
        Nuc_vid(:,:,i) = nuimage;
    end

    Islet_vid = rot90(Islet_vid, 2);%the image is reflected rotate it back
    Nuc_vid = rot90(Nuc_vid, 2);
    ca_im = rescale(mean(Islet_vid(:,:,timetouse(kt)-10:timetouse(kt)+10),3),0,1);%-double(min(Islet_vid,[],3)))./(double(max(Islet_vid,[],3))-double(min(Islet_vid,[],3)));%mean2([squeeze(mean(Islet_vid,2)),squeeze(mean(Islet_vid,1))]));
    
    


    try 
       load(strrep(strjoin([savepath filename(kt) '_' masktype '.mat']), ' ', ''))
    catch
        %% Load nucleus location 
        nucloc = readtable(strrep(strjoin([datapath csvname(kt)]),'/ ','/')); %import nucleus location csv
        cells_at_time = find(table2array(nucloc(:,7))==timetouse(kt));
        X = (table2array(nucloc(cells_at_time,1)));
        Y = table2array(nucloc(cells_at_time,2));

        nuimage = Nuc_vid(:,:,timetouse(kt));
        loc = [X,Y]; 

        trainingcells = randi([1,length(X)], 1,perislet); %select 'perislet' number of training cells at random with generated seed


        imnew = insertMarker(nuimage, loc);
        figure, imshow(imnew)

        cells = [1:length(X)]; %array of cells
        numcells = length(X);
        %% Grab nucleus location and cell outline
        %create a color plot of nucleus and calcium
        caim_nuc = insertMarker(ca_im, loc(:,:));
        fig = figure('Name','Draw ROIs Around Cells of Interest');
        fig.Position = [1293 366 1076 918];
        fpos = get(fig,'Position')
        axOffset = (fpos(3:4)-[size(caim_nuc,2) size(caim_nuc,1)])/2;
        ha = axes('Parent',fig,'Units','pixels',...
                    'Position',[axOffset size(caim_nuc,2) size(caim_nuc,1)]);
        imshow(caim_nuc, 'Parent',ha);
        CellMask = double(zeros(size(ca_im)));
        for i = trainingcells
            imshow(caim_nuc)

            NucLoc(i,:) = loc(i,:);
            caim_onenuc = insertMarker(ca_im, NucLoc(i,:));
            imshow(caim_onenuc)
            title('Draw Around Cell')
          
            ROIMask = imfreehand(); %User draws region around cell
            ROIMask = createMask(ROIMask); %Mask is created from drawn region
            CellMask = CellMask + ROIMask.*i; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
            caim_nuc = insertMarker(caim_nuc, NucLoc(i,:), 'color','r');
        end
        imshow(caim_nuc)
        NucLoc = fliplr(NucLoc);

    end
    
   %Show image with cell masks 
    cells_outline = imfuse(caim_nuc, CellMask);
    cells_w_labels = figure;
    imshow(cells_outline)
    title(['Islet ' filename(kt)])

    for c = 1:perislet
        text(loc(trainingcells(c),1), loc(trainingcells(c),2), num2str(c), 'Color', 'w'); % Labels cells in the image with their respective region number
    end

    saveas(cells_w_labels, (strrep(strjoin([savepath '/Figures/Masks_ ' filename(kt) '_' masktype '.fig']), ' ', ''))); % Saves connection map
    saveas(cells_w_labels, (strrep(strjoin([savepath '/Figures/Masks_ ' filename(kt) '_' masktype '.png']), ' ', ''))); % Saves connection map
    

    %Extract Calcium:
    for i = 1:perislet
        TCMask = CellMask;
        TCMask(CellMask ~= trainingcells(i)) = 0;
        MaskedIMGstack = double(Islet_vid).*logical(TCMask);
        [r,c] = (find(MaskedIMGstack(:,:,10))); %find mask at time == 10 just in case there is something wrong with the first frame
        for k = 1:length(r)
            Ca(k,:) = MaskedIMGstack(r(k),c(k),:);
        end

        Calcium(:,i) = mean(Ca);
        clear Ca
    end
    
    %Calculate Correlation & Network Info:

    %First do 5 links
    Opts.Method = 'Degree';
    Opts.avDeg = 5;
    Opts.figs = 0;
    Thr = findoptRth(Calcium, Opts)
    [N, Adj, kpercent, histArrayPercShort,pval,Rij,s] = NetworkAnalysis(Calcium, Thr, Opts,0)%ii, mm, phase, figs)
    

    %ST analysis:

    %Calculate Correlatino & Network Info:


    save(strrep(strjoin(['/Volumes/Briggs_10TB/Merrin/Confocal/' filename(kt) 'NormalizeAnalysis.mat']), ' ', ''), 'CellMask','FinalCordata','NucLoc', 'trainingcells')


end
