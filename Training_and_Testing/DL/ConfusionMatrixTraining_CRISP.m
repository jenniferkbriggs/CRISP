%Load masks from file:
clearvars -except confMatper
close all

files = dir('/Volumes/Briggs_10TB/CRISPdata/Analysis/*.csv')

addpath('~/Documents/GitHub/CRISP/Run_CRISP/')
savepath = '/Users/brigjenn/Documents/GitHub/ST_Analysis/Data/'

%load reference mask:
load(([savepath 'two_Good.mat']))
reference=CellMask;
%load DL mask:
load(([savepath 'two_Medium.mat']))

Opts.fig =0;
Opts.Thr = 'st';
Opts.st_thr=0.5;
tic
CellMask = Mask_refinement(Islet_vid, CellMask, Opts)
toc

%now i have to go through and find indexes close to the cell mask.
%NucLockeep


caim_mask = insertMarker(CellMask, ((NucLoc_keep)),"circle","Color","red");
figure2=figure
imshow(caim_mask)


    %negative = remove, positive = keep
RefernenceB=reference>0;
CellMaskB=CellMask>0;

refVec=RefernenceB(:);
cellVec = CellMaskB(:);

confMat=confusionmat(refVec,cellVec);
totalPix=sum(sum(confMat));
confMatper(:,:,5)=confMat./totalPix;

 
%% - once you are done with all of them:
if length(find(confMatper==0))==0
    
    confusionMatrix = mean(confMatper,3);
    confusionMatrix_std=std(confMatper,[],3)
    figure;
    set(gcf,'Position',[817,788,743,519]);
    set(gcf,"Units",'pixels')
    imagesc(confusionMatrix);
    
    % Set custom colormap: white for negatives, green for positives
    % customCMap = [1 1 1; 79./255 121./255 66./255]; % White for negatives, green for positives
    % colormap(customCMap);
    
    % Add text annotations with mean ± std
    textStrings = strcat(num2str(confusionMatrix(:),'%.3f'), '\pm', num2str(confusionMatrix_std(:),'%.3f'));
    textStrings = strtrim(cellstr(textStrings));
    [x, y] = meshgrid(1:2);
    hStrings = text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center', 'FontSize', 25);
    set(gca,'ColorScale','log')
    colormap(flipud(autumn))

    % Set axis labels and other properties
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Negative', 'Positive'}, ...
             'YTick', 1:2, 'YTickLabel', {'Negative', 'Positive'}, ...
             'TickLength', [0 0]);
    xlabel('Predicted Class');
    ylabel('Actual Class');
    set(gca, 'FontSize', 30);
    set(gcf, 'Color', 'white');


saveas(gcf, '/Volumes/Briggs_10TB/CRISPdata/Analysis/TrainingConfusion_pixel_cellpose.fig')
saveas(gcf, '/Volumes/Briggs_10TB/CRISPdata/Analysis/TrainingConfusion_pixel_cellpose.png')
save('/Volumes/Briggs_10TB/CRISPdata/Analysis/TrainingConfusion_pixel_cellpose.mat', 'confMatper')





end