%%Estimate radius of all cells in the nuclei
%close all
clear all
clc

opts.score_thr = 0.9
opts.th_pix =  0.8183

addpath('~/Documents/GitHub/UniversalCode/');

%% Set up files to analyze: chose islet to analyze
filename = ["three","sample","five", "two","one"];
csvname = ["H2BmCherry Ucn3GCaMP-3_Detailed.csv", "H2BmCherry Ucn3GCaMP sample_Detailed.csv", "H2BmCherry Ucn3GCaMP-5_Detailed.csv", "H2BmCherry Ucn3GCaMP-2_Detailed.csv", "H2BmCherry Ucn3GCaMP-1_Detailed.csv"]
datapath = ['/Volumes/Briggs_10TB/Merrin/Confocal/'] 

maxradius = 40;

circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

%initialize which time index to use nuclear locations
%timetouse = 243%three: 217 - somewhat wiggly. %sample: 194; %five: 243%two: 363; %one: 193
timetouse = [217, 194, 243, 363, 193];

for kt = 1:length(filename)
    isletfig = figure

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
    ca_im = rescale(mean(Islet_vid(:,:,timetouse(kt)-10:timetouse(kt)+10),3),0,1);

    imshow(ca_im)
   %% Load nucleus location 
        nucloc = readtable(strrep(strjoin([datapath csvname(kt)]),'/ ','/')); %import nucleus location csv
        cells_at_time = find(table2array(nucloc(:,7))==timetouse(kt));
        X = (table2array(nucloc(cells_at_time,1)));
        Y = table2array(nucloc(cells_at_time,2));

        nuimage = Nuc_vid(:,:,timetouse(kt));
        loc = [X,Y]; 

        imnew = insertMarker(ca_im, loc(61,:));
        figure, imshow(imnew)
        
        NucLoc = fliplr(loc);

    %Because of that rotation, my nucleus location is rotated: flip back
    imAv = mean(Nuc_vid, 3);
    imAv = imAv/max(imAv(:));
    for i = 61:length(NucLoc)
        opts.ca = ca_im;
        opts.gif = 0;
        opts.normalized = 1;
        opts.title = [filename(kt) 'Islet' num2str(i)]
        
        %HERE IS WHERE CRISP IS RUN --- 

       [Correlation, radius, pixelsx, pixelsy] = CRISP_annulus_corr(images,NucLoc,opts)
        est_radius(kt,i) = radius; 
end 
end

save([datapath 'AllCellResults.mat'], 'est_radius')
