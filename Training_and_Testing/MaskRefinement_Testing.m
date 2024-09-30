%%Here we mask and test testing data for cell mask refinement:
close all
clear all
clc

addpath('~/Documents/GitHub/UniversalCode/');
addpath('~/Documents/GitHub/Islet_Heterogeneity/')
%chose islets to analyze
pathname = '/Volumes/Briggs_10TB/NetworkPaper/22nd Aug/22nd Aug/'

ending = '.lsm';
savename = '/PlaywithSTD.mat'
numad = 1;
zstacks = 1;
zz = 1;
savepath = '/Users/brigjenn/OneDrive - The University of Colorado Denver/Anschutz/Islet/STDanalysis/'

cachannel = 1;
howmanychannel = 1;

filename = dir([pathname, 'Cx36mm*']);
savepath = ['~/Documents/GitHub/ST_Analysis/Data/']



%initialize which time index to use nuclear locations
%timetouse = 243%three: 217 - somewhat wiggly. %sample: 194; %five: 243%two: 363; %one: 193
timetouse = [217, 194, 243, 363, 193];
%Number of cells per islet for the training case
perislet = 10;

%set seed

for kt = [1:3]
for mm = [1:3]
    %load:
    R = bfopen([pathname filename(kt).name]); % Uses bfopen program to open .czi/.lsm image files
    try
        save([savepath savename],  '-V7.3')
    catch
        mkdir([savepath])
        save([savepath savename],  '-V7.3')
    end


%set random number generator
rng(kt*2)



tic
pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end

try
    for i=1:length(pics)
        T(i)=R{4}.getPlaneDeltaT(0, i-1).value;
    end
catch
    T=0:0.5:pn*0.5;
end
T = double(T);
T = T(cachannel:howmanychannel:end);
T = T(1:zstacks:end);
%T = T(3:3:end);

    st = 1;

    ed=length(T);


T = T(st:ed);

images=double(IMG); % converts images to double precision
images = images(:,:,cachannel:howmanychannel:end);

RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

clear pics R IMG;
% clear omeMeta;

output_dir = savepath;
toc

%% DECLARING IMAGE PROPERTIES
tic
images = images(:,:,zz:zstacks:end);
images = images(:,:,st:ed-1);
sx=size(images,1);
sy=size(images,2);
sz=length(T)-1;
for i=1:sz
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
end

%Makes all frames 1fps
% fps =Vidinfo(g).fps(ff);
fps = 1;
 ct =1
for i = 1:fps:sz
    images(:,:,ct) = images(:,:,i);
    ct = ct+1;
end
images = images(:,:,1:ct-1);
%ImAv = sum(images,3); %compresses all frames into single array of intensities
ImAv = mean(images,3); 
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image
Islet_vid = images;
RGB2 = hsv2rgb(HSV); %converts to rgb image

RGB= im2gray(RGB2);
RGB = imadjust(RGB);
OGFig = figure(1);
caim_nuc = RGB;
ca_im = RGB;
imshow(RGB)


    masktypes = {'Bad', 'Medium', 'Good'}
    masktype = masktypes{mm}
    %% Import Images %%
    %import calcium
   

    try 
       load((([savepath filename(kt).name '_' masktype '.mat'])),...
           'CellMask', 'RGB', 'caim_nuc', 'trainingcells', 'Islet_Vid', 'NucLoc')


    catch
        %% Load nucleus location 

        cells = [1:perislet]; %array of cells
        numcells = perislet;
        %% Grab nucleus location and cell outline
        %create a color plot of nucleus and calcium

        CellMask = double(zeros(size(caim_nuc)));

        if mm == 1
            for i = 1:perislet
            %for first mask type, choose nuclei 
            % Create the figure and set its properties
            fig = figure;
            fig.Position = [1607 527 737 658]; % Set the position and size of the figure
            fig.Units = 'pixels'; % Set units to pixels
            
            % Display the image in the figure
            imshow(caim_nuc);
            
            % Get the current axes
            ax = gca;
            
            % Set the axes properties to make the image larger
            ax.Position = [0 0 1 1]; % Make the axes fill the entire figure
            axis on; % Optionally turn on axis labels if needed
            
            % Optionally, you can adjust the image aspect ratio if desired
            axis image; % Ensure the aspect ratio of the image is preserved
            fig.Position = [1607 527 737 658]; % Set the position and size of the figure
            fig.Units = 'pixels'; % Set units to pixels


           
                title('Select Nucleus')
                training_loc = ginput(1);
                trainingcells(i) = i;
                NucLoc(i,:) = training_loc;
                caim_onenuc = insertMarker(ca_im, NucLoc(i,:));

                           fig = figure;
            fig.Position = [1607 527 737 658]; % Set the position and size of the figure
            fig.Units = 'pixels'; % Set units to pixels
            
            % Display the image in the figure
            imshow(caim_onenuc);
            
            % Get the current axes
            ax = gca;
            
            % Set the axes properties to make the image larger
            ax.Position = [0 0 1 1]; % Make the axes fill the entire figure
            axis on; % Optionally turn on axis labels if needed
            
            % Optionally, you can adjust the image aspect ratio if desired
            axis image; % Ensure the aspect ratio of the image is preserved
            fig.Position = [1607 527 737 658]; % Set the position and size of the figure
            fig.Units = 'pixels'; % Set units to pixels

                title('Zoom then press continue')
                keyboard
                title('Draw around cell')
                ROIMask = imfreehand(); %User draws region around cell
                ROIMask = createMask(ROIMask); %Mask is created from drawn region
                CellMask = CellMask + ROIMask.*i; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
                CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
                caim_nuc = insertMarker(caim_nuc, NucLoc(i,:), 'color','r');
                disp(num2str(i))
                            CellMask(CellMask > perislet) = 0;

            end
        else
            for i = 1:perislet
                disp(masktypes{mm})
                %image flips but I save mm = 2 as fliped so mm = 3 doesn't
                %need flipped

                caim_onenuc = insertMarker(ca_im, ((NucLoc(i,:))));
                fig = figure;
                fig.Position = [1607 527 737 658]; % Set the position and size of the figure
                fig.Units = 'pixels'; % Set units to pixels
                
                % Display the image in the figure
                imshow(caim_onenuc);
                
                % Get the current axes
                ax = gca;
                
                % Set the axes properties to make the image larger
                ax.Position = [0 0 1 1]; % Make the axes fill the entire figure
                axis on; % Optionally turn on axis labels if needed
                
                % Optionally, you can adjust the image aspect ratio if desired
                axis image; % Ensure the aspect ratio of the image is preserved
                fig.Position = [1607 527 737 658]; % Set the position and size of the figure
                fig.Units = 'pixels'; % Set units to pixels


                title('Zoom then press continue')
                keyboard
                title('Draw around cell')
                ROIMask = imfreehand(); %User draws region around cell
                ROIMask = createMask(ROIMask); %Mask is created from drawn region
                CellMask = CellMask + ROIMask.*i; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
                CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
                caim_nuc = insertMarker(caim_nuc, NucLoc(i,:), 'color','r');
                disp(num2str(i))
            end
            %Get rid of any other overlap:
            CellMask(CellMask > perislet) = 0;
        end
        imshow(caim_nuc)
        close all
        save([savepath filename(kt).name '_' masktype '.mat'])


    end

  
   %Show image with cell masks 
    cells_outline = imfuse(caim_nuc, CellMask);
    cells_w_labels = figure;
    imshow(cells_outline)
    title(['Islet ' filename(kt).name])

    for c = 1:perislet
        text(NucLoc(trainingcells(c),1), NucLoc(trainingcells(c),2), num2str(c), 'Color', 'w'); % Labels cells in the image with their respective region number
    end

    saveas(cells_w_labels, ((([savepath '/Figures/Masks_ ' filename(kt).name '_' masktype '.fig'])))); % Saves connection map
    saveas(cells_w_labels, ((([savepath '/Figures/Masks_ ' filename(kt).name '_' masktype '.png'])))); % Saves connection map
    
    if 1
    Calcium = [];
    %Extract Calcium:
    for i = 1:perislet
        TCMask = CellMask;
        TCMask(CellMask ~= (i)) = 0;
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
    
   % load([savepath 'RefineMasks.mat'],'output'); % Saves connection map


    
    output.(filename(kt).name(1:end-4)).(masktype).corr = mean(nonzeros(triu(Rij,1)));
    output.(filename(kt).name(1:end-4)).(masktype).Rij = Rij;
    output.(filename(kt).name(1:end-4)).(masktype).N = N;
    output.(filename(kt).name(1:end-4)).(masktype).Calcium = Calcium;
    output.(filename(kt).name(1:end-4)).(masktype).Adj = Adj;
    output.(filename(kt).name(1:end-4)).(masktype).Thr = Thr;
    output.(filename(kt).name(1:end-4)).(masktype).CellMask_st = CellMask;

    %ST analysis:
    %Opts: 
    Opts.fig = 1;
    Opts.st_thr =0.25;
    Opts.Thr = 'st'
        CellMask = Mask_refinement(Islet_vid, CellMask, Opts);
        saveas(gcf, ((([savepath '/Figures/MasksRefinedSTAnalysis_ ' filename(kt).name '_' masktype '.fig'])))); % Saves connection map
        saveas(gcf, ((([savepath '/Figures/MasksRefinedSTAnalysis_ ' filename(kt).name '_' masktype '.png'])))); % Saves connection map
        
    %Extract Calcium:
    for i = 1:perislet
        TCMask = CellMask;
        TCMask(CellMask ~= (i)) = 0;
        MaskedIMGstack = double(Islet_vid).*logical(TCMask);
        [r,c] = (find(MaskedIMGstack(:,:,10))); %find mask at time == 10 just in case there is something wrong with the first frame
        for k = 1:length(r)
            Ca(k,:) = MaskedIMGstack(r(k),c(k),:);
        end

        try
        Calcium_2(:,i) = mean(Ca);
        end
        clear Ca
    end
    
    %Calculate Correlatin & Network Info:
    %Thr = findoptRth(Calcium_2, Opts) - use threshold from before:
    [N, Adj, kpercent, histArrayPercShort,pval,Rij,s] = NetworkAnalysis(Calcium_2, Thr, Opts,0)%ii, mm, phase, figs)
    
    output.(filename(kt).name(1:end-4)).(masktype).corr_st = mean(nonzeros(triu(Rij,1)));
    output.(filename(kt).name(1:end-4)).(masktype).Rij_st = Rij;
    output.(filename(kt).name(1:end-4)).(masktype).N_st = N;
    output.(filename(kt).name(1:end-4)).(masktype).Calcium_st = Calcium_2;
    output.(filename(kt).name(1:end-4)).(masktype).Adj_st = Adj;
    output.(filename(kt).name(1:end-4)).(masktype).Thr_st = Thr;
    output.(filename(kt).name(1:end-4)).(masktype).CellMask_st = CellMask;

    save([savepath 'RefineMasks.mat'],'output'); % Saves connection map


    clear Calcium2
    end
end
end
close all


%% -- 
filename = fieldnames(output)

% Analysis: 

for i = 1:3
file = string(filename(i));
Bad(i) = output.(file).Bad.corr;
Bad_s(i) = output.(file).Bad.corr_st;
Med(i) = output.(file).Medium.corr;
Med_s(i) = output.(file).Medium.corr_st;
Good(i) = output.(file).Good.corr;
Good_s(i) = output.(file).Good.corr_st;
end
writematrix(([Bad,Bad_s; Med, Med_s]), [savepath, 'Corr.csv'])

for i = 1:3
%difference between Rij: 
file = string(filename(i));
output.(file).Bad.Rij(logical(eye(size(output.(file).Bad.Rij)))) = NaN;
output.(file).Bad.Rij_st(logical(eye(size(output.(file).Bad.Rij)))) = NaN;
output.(file).Good.Rij_st(logical(eye(size(output.(file).Bad.Rij)))) = NaN;
output.(file).Good.Rij(logical(eye(size(output.(file).Bad.Rij)))) = NaN;
output.(file).Medium.Rij(logical(eye(size(output.(file).Bad.Rij)))) = NaN;
output.(file).Medium.Rij_st(logical(eye(size(output.(file).Bad.Rij)))) = NaN;
end

for i = 1:3
file = string(filename(i));

Bad(i) = mean(mean(output.(file).Bad.Rij-output.(file).Good.Rij, 'omitnan'), 'omitnan');
Bad_s(i) = mean(mean(output.(file).Bad.Rij_st-output.(file).Good.Rij, 'omitnan'), 'omitnan');
Med(i) = mean(mean(output.(file).Medium.Rij- output.(file).Good.Rij, 'omitnan'), 'omitnan');
Med_s(i) = mean(mean(output.(file).Medium.Rij_st-output.(file).Good.Rij, 'omitnan'), 'omitnan');
end
writematrix(([Bad,Bad_s; Med, Med_s]), [savepath, 'Rij.csv'])




for i = 1:3
file = string(filename(i));
Bad(i) = max(output.(file).Bad.N);
Bad_s(i) = max(output.(file).Bad.N_st);
Med(i) = max(output.(file).Medium.N);
Med_s(i) = max(output.(file).Medium.N_st);
Good(i) = max(output.(file).Good.N);
Good_s(i) = max(output.(file).Good.N_st);
end
writematrix(([Bad', Med', Good', Bad_s', Med_s', Good_s']')', [savepath, 'Max_N.csv'])


for i = 1:3
file = string(filename(i));
Bad(i) = mean(output.(file).Bad.N);
Bad_s(i) = mean(output.(file).Bad.N_st);
Med(i) = mean(output.(file).Medium.N);
Med_s(i) = mean(output.(file).Medium.N_st);
Good(i) = mean(output.(file).Good.N);
Good_s(i) = mean(output.(file).Good.N_st);
end
writematrix(([Bad', Med', Good', Bad_s', Med_s', Good_s']')', [savepath, 'Mean_Degree.csv'])


for i = 1:3
file = string(filename(i));
figure,
histogram(output.(file).Bad.N, 'numbins',5);
hold on, histogram(output.(file).Good.N, 'numbins',5);
legend('Bad','Good')

figure,
histogram(output.(file).Bad.N_st, 'numbins',5);
hold on, histogram(output.(file).Good.N, 'numbins',5);
legend('Bad_st','Good')

figure,
histogram(output.(file).Good.N_st, 'numbins',5);
hold on, histogram(output.(file).Good.N, 'numbins',5);
legend('Good_st','Good')

figure,
histogram(output.(file).Medium.N, 'numbins',5);
hold on, histogram(output.(file).Good.N, 'numbins',5);
legend('Medium','Good')

figure,
histogram(output.(file).Medium.N_st, 'numbins',5);
hold on, histogram(output.(file).Good.N, 'numbins',5);
legend('Medium_st','Good')


end



