%Load masks from file:
clearvars -except confMatper
close all

files = dir('/Volumes/Briggs_10TB/CRISPdata/Analysis/*.csv')


savepath = '/Users/brigjenn/Documents/GitHub/ST_Analysis/Data/'
% %load reference mask:
% load(([savepath files(10).name(2:end-9) '.lsm_Good.mat']))
% reference=CellMask;
% 
% %load DL mask:
% CellMask = readmatrix(['/Volumes/Briggs_10TB/CRISPdata/Analysis/' files(10).name]);


%load reference mask:
load(([savepath 'five_Good.mat']))
reference=CellMask;
ll = 12;
%load DL mask:
CellMask = readmatrix(['/Volumes/Briggs_10TB/CRISPdata/Analysis/' files(ll).name]);
%exports in python so first index is not real
CellMask(1,:)= [];


%now i have to go through and find indexes close to the cell mask.
%NucLockeep


caim_mask = insertMarker(CellMask, ((NucLoc_keep)),"circle","Color","red");
figure2=figure
imshow(caim_mask)

%fin
    for c = 1:max(max(CellMask))
        [x,y] = find(CellMask ==c)
        text(mean(y),mean(x), num2str(c)); % Labels cells in the image with their respective region number
    end



ct=1;
k=1
while k ==1
[x, y]=ginput(1);
maskvalue(ct) = CellMask(round(y),round(x))
ct = ct+1;
if mod(ct,3)==0
k = input("keep going?")
end
end

%adjust as needed
keyboard

%remove Pixels in mask:
for i =unique(maskvalue)
    [x,y] =find(CellMask == i);
    for j= 1:length(x)
        CellMask(x(j),y(j)) = 0;
    end
end

%save:



    %negative = remove, positive = keep
RefernenceB=reference>0;
CellMaskB=CellMask>0;

refVec=RefernenceB(:);
cellVec = CellMaskB(:);

confMat=confusionmat(refVec,cellVec);
totalPix=sum(sum(confMat));
confMatper(:,:,5)=confMat./totalPix;
save(['/Volumes/Briggs_10TB/CRISPdata/Analysis/' files(ll).name 'MaskWithReferenceCells.mat']);

 
%% - once you are done with all of them:
if length(find(confMatper==0))==0
    
    confusionMatrix = mean(confMatper,3);
    confusionMatrix_std=std(confMatper,[],3)
    figure;
    set(gcf,'Position',[817,788,743,519]);
    set(gcf,"Units",'px')
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
    colormap(flipud(summer))

    % Set axis labels and other properties
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Negative', 'Positive'}, ...
             'YTick', 1:2, 'YTickLabel', {'Negative', 'Positive'}, ...
             'TickLength', [0 0]);
    xlabel('Predicted Class');
    ylabel('Actual Class');
    set(gca, 'FontSize', 30);
    set(gcf, 'Color', 'white');

saveas(gcf, '/Volumes/Briggs_10TB/CRISPdata/Analysis/TrainingConfusion_pixel.fig')
saveas(gcf, '/Volumes/Briggs_10TB/CRISPdata/Analysis/TrainingConfusion_pixel.png')


end

% Accuracy: 
for i = 1:5
    %should already be normalized (per means percent)
    if (confMatper(1,1,i) + confMatper(2,2,i)+confMatper(1,2,i) + confMatper(2,1,i)) ~= 1
        disp('Error')
        keyboard
    end
    acc(i) = (confMatper(1,1,i) + confMatper(2,2,i));
    sens(i) = confMatper(1,1,i)./(confMatper(1,1,i)+confMatper(2,1,i));
    spec(i) = confMatper(2,2,i)./(confMatper(2,2,i)+confMatper(1,2,i));
end


