%%Here we mask and training data set for cell mask refinement:
close all
clear all
clc

addpath('~/Documents/GitHub/UniversalCode/');
addpath('~/Documents/GitHub/Islet_Heterogeneity/')
%chose islets to analyzeclose all

cachannel = 1;
nuchannel = 2;
howmanychannel = 1;

savepath = '/Users/brigjenn/Documents/GitHub/ST_Analysis/Data/'

%chose islet to analyze
filename = ["three","five", "two"]%,"sample","one"];
csvname = ["H2BmCherry Ucn3GCaMP-3_Detailed.csv", "H2BmCherry Ucn3GCaMP-5_Detailed.csv", "H2BmCherry Ucn3GCaMP-2_Detailed.csv"]%,...
    %"H2BmCherry Ucn3GCaMP sample_Detailed.csv",  "H2BmCherry Ucn3GCaMP-1_Detailed.csv"]
datapath = ['/Volumes/Briggs_10TB/Merrin/Confocal/'] 


%initialize which time index to use nuclear locations
%timetouse = 243%three: 217 - somewhat wiggly. %sample: 194; %five: 243%two: 363; %one: 193
timetouse = [217, 194, 243, 363, 193];
%Number of cells per islet for the training case
perislet = 15;

%set seed

for kt = 1:3

%set random number generator
rng(kt*2)


for mm = [1:3]


ca_files = dir(strrep(strjoin([datapath filename(kt) '/' '*C1*.tif']),' ',''));
nuc_files = dir(strrep(strjoin([datapath filename(kt) '/' '*C2*.tif']),' ',''));

    masktypes = {'Bad', 'Medium', 'Good'}
    masktype = masktypes{mm}
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
       load(strrep(strjoin([savepath filename(kt) '_' masktype '.mat']), ' ', ''),'loc',...
           'CellMask', 'Islet_vid', 'caim_nuc', 'trainingcells')


    catch
        %% Load nucleus location 
        nucloc = readtable(strrep(strjoin([datapath csvname(kt)]),'/ ','/')); %import nucleus location csv
        cells_at_time = find(table2array(nucloc(:,7))==timetouse(kt));
        X = (table2array(nucloc(cells_at_time,1)));
        Y = table2array(nucloc(cells_at_time,2));

        nuimage = Nuc_vid(:,:,timetouse(kt));
        loc = [X,Y]; 


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

        if mm == 1
            %for first mask type, choose nuclei 
            for i = 1:perislet
                imshow(caim_nuc)
                title('Zoom and click continue')
                keyboard
                title('Select Nucleus')
                training_loc = ginput(1);
                NucLoc_Ithink = find(sum(abs(loc - training_loc),2)<5);
                if length(NucLoc_Ithink)~= 1
                    keyboard
                end
                trainingcells(i) = NucLoc_Ithink;
                NucLoc(i,:) = loc(NucLoc_Ithink,:);
                caim_onenuc = insertMarker(ca_im, NucLoc(i,:));
                imshow(caim_onenuc)
                title('Zoom then press continue')
                keyboard
                title('Draw around cell')
                ROIMask = imfreehand(); %User draws region around cell
                ROIMask = createMask(ROIMask); %Mask is created from drawn region
                CellMask = CellMask + ROIMask.*i; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
                CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
                caim_nuc = insertMarker(caim_nuc, NucLoc(i,:), 'color','r');
                disp(num2str(i))
                            CellMask(CellMask > 15) = 0;

            end
        else
            for i = 1:perislet
                disp(masktypes{mm})
                %image flips but I save mm = 2 as fliped so mm = 3 doesn't
                %need flipped

                caim_onenuc = insertMarker(ca_im, ((NucLoc(i,:))));

                imshow(caim_onenuc)
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
            CellMask(CellMask > 15) = 0;
        end
        imshow(caim_nuc)
        save(strrep(strjoin([savepath filename(kt) '_' masktype '.mat']), ' ', ''))


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

    %load([savepath 'RefineMasks.mat'],'output'); % Saves connection map

    output.(filename(kt)).(masktype).corr = mean(nonzeros(triu(Rij,1)));
    output.(filename(kt)).(masktype).N = N;
    output.(filename(kt)).(masktype).Calcium = Calcium;
    output.(filename(kt)).(masktype).Adj = Adj;
    output.(filename(kt)).(masktype).Thr = Thr;
    output.(filename(kt)).(masktype).CellMask_st = CellMask;
    output.(filename(kt)).(masktype).Rij = Rij;
    output.(filename(kt)).(masktype).averagephase = averagephase;
    output.(filename(kt)).(masktype).sorted_highphase = sorted_highphase;

    %ST analysis:
        %Opts: 
    Opts.fig = 0;
    Opts.st_thr =0.25;
    Opts.Thr = 'st'
        CellMask = Mask_refinement(Islet_vid, CellMask, Opts);
        saveas(gcf, (strrep(strjoin([savepath '/Figures/MasksRefinedSTAnalysis_ ' filename(kt) '_' masktype '.fig']), ' ', ''))); % Saves connection map
        saveas(gcf, (strrep(strjoin([savepath '/Figures/MasksRefinedSTAnalysis_ ' filename(kt) '_' masktype '.png']), ' ', ''))); % Saves connection map
        
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
    [averagephase, sorted_highphase]  = RunPhaseAnalysis_allsecondphase(Calcium_2);

    output.(filename(kt)).(masktype).corr_st = mean(nonzeros(triu(Rij,1)));
    output.(filename(kt)).(masktype).N_st = N;
    output.(filename(kt)).(masktype).Calcium_st = Calcium_2;
    output.(filename(kt)).(masktype).Adj_st = Adj;
    output.(filename(kt)).(masktype).Thr_st = Thr;
    output.(filename(kt)).(masktype).CellMask_st = CellMask;
        output.(filename(kt)).(masktype).averagephase_st = averagephase;
        output.(filename(kt)).(masktype).sorted_highphase_st = sorted_highphase;


        output.(filename(kt)).(masktype).Rij_st = Rij;

    save([savepath 'RefineMasks.mat'],'output'); % Saves connection map


    clear Calcium2 
    end
end
end
close all

% Analysis: 

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

Bad2(i,:) = mean(output.(file).Bad.Rij-output.(file).Good.Rij, 'omitnan');
Bad_s2(i,:) = (mean(output.(file).Bad.Rij_st-output.(file).Good.Rij, 'omitnan'));
Med2(i,:) = (mean(output.(file).Medium.Rij- output.(file).Good.Rij, 'omitnan'));
Med_s2(i,:) = (mean(output.(file).Medium.Rij_st-output.(file).Good.Rij, 'omitnan'));
end
writematrix(([reshape(Bad2, [],1),reshape(Bad_s2, [],1)]), [savepath, 'Rij_manybad.csv'])
writematrix(([reshape(Med2,[],1), reshape(Med_s2,[],1)]), [savepath, 'Rij_manymed.csv'])





for i = 1:3
file = string(filename(i));

Bad(i) = mean(mean(abs(output.(file).Bad.averagephase-output.(file).Good.averagephase), 'omitnan'), 'omitnan');
Bad_s(i) = mean(mean(abs(output.(file).Bad.averagephase_st-output.(file).Good.averagephase), 'omitnan'), 'omitnan');
Med(i) = mean(mean(abs(output.(file).Medium.averagephase- output.(file).Good.averagephase), 'omitnan'), 'omitnan');
Med_s(i) = mean(mean(abs(output.(file).Medium.averagephase_st-output.(file).Good.averagephase), 'omitnan'), 'omitnan');
end
writematrix(([Bad,Bad_s; Med, Med_s]), [savepath, 'PhaseDiff.csv'])


for i = 1:3
file = filename(i);
Bad(i) = output.(file).Bad.corr;
Bad_s(i) = output.(file).Bad.corr_st;
Med(i) = output.(file).Medium.corr;
Med_s(i) = output.(file).Medium.corr_st;
Good(i) = output.(file).Good.corr;
Good_s(i) = output.(file).Good.corr_st;
end
writematrix([Bad'-Good',Bad_s'-Good'], [savepath, 'Corr_Bad.csv'])

writematrix(([Med'-Good', Med_s'-Good']), [savepath, 'Corr_Med.csv'])

for i = 1:3
file = filename(i);
Bad(i) = max(output.(file).Bad.N);
Bad_s(i) = max(output.(file).Bad.N_st);
Med(i) = max(output.(file).Medium.N);
Med_s(i) = max(output.(file).Medium.N_st);
Good(i) = max(output.(file).Good.N);
Good_s(i) = max(output.(file).Good.N_st);
end
writematrix(([Bad',Bad_s', Med', Med_s']-Good')', [savepath, 'Max_N.csv'])


for i = 1:4
file = filename(i);
Bad(i) = mean(output.(file).Bad.N);
Bad_s(i) = mean(output.(file).Bad.N_st);
Med(i) = mean(output.(file).Medium.N);
Med_s(i) = mean(output.(file).Medium.N_st);
Good(i) = mean(output.(file).Good.N);
Good_s(i) = mean(output.(file).Good.N_st);
end
writematrix(([Bad', Med', Good', Bad_s', Med_s', Good_s']./Good')', [savepath, 'Mean_Degree.csv'])


for i = 1:4
file = filename(i);
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



