%% Veroni tesselation: 

close all
clearvars -except confMatper
clc

%Find the radial location of training
th_pix = 0.8183; %threshold for pixels
th_rad = 0.83; %score for radius. 

indx = 5
%load:
savepath = '/Users/brigjenn/Documents/GitHub/ST_Analysis/Data/'
%load data
title = 'five'
load(([savepath 'five_Good.mat']))
reference = CellMask;
%voronoi: 
[v,c] = voronoin(loc);
%i need to write a code that converts this into a cell mask type file:
CellMask = zeros(size(CellMask));
badindx= find(isinf(v))
for ci = 1:length(c)
    %retrieve edge index:
    edges = c{ci};
    %if edge is Inf - then it reaches the edge of the screen. just omit: 
    for j = 1:length(badindx)
        edges(find(edges == badindx(j))) = [];
    end
    incell = poly2mask(v(edges,1), v(edges,2), size(CellMask,1), size(CellMask,2));
    CellMask=CellMask + incell.*ci;
end

rgbim = label2rgb(CellMask, 'jet','k','shuffle'); figure, imshow(rgbim)
hold on, plot(loc(:,1), loc(:,2), 'o', 'Color','k', 'MarkerSize',4)

saveas(gcf, ['/Volumes/Briggs_10TB/CRISPdata/Analysis/Voronoi_' title '.png'])

%% masks indexes should already overlap. 
caim_mask = insertMarker(CellMask, ((NucLoc_keep)),"circle","Color","red");
figure2=figure
imshow(caim_mask)


    for c = 1:max(max(CellMask))
        [x,y] = find(CellMask ==c)
        text(mean(y),mean(x), num2str(c)); % Labels cells in the image with their respective region number
    end


keyboard
ct=1;
k=1
while k ==1
[x, y]=ginput(1);
maskvalue(ct) = CellMask(round(y),round(x))
if mod(ct,15)==0
k = input("keep going?")
end
ct = ct+1;

end

%adjust as needed
keyboard

CellMaskNew = zeros(size(CellMask,1), size(CellMask,2))
%remove Pixels in mask:
for i =unique(maskvalue)
    [x,y] =find(CellMask == i);
    for j= 1:length(x)
        CellMask_New(x(j),y(j)) = i;
    end
end

v_outline = imfuse(ca_im, CellMask_New);
figure, imshow(v_outline)
saveas(gcf, ['/Volumes/Briggs_10TB/CRISPdata/Analysis/Voronoi_onlyselectedcells' title '.png'])


%%compute confusion mathrix:



    %negative = remove, positive = keep
RefernenceB=reference>0;
CellMaskB=CellMask>0;

refVec=RefernenceB(:);
cellVec = CellMaskB(:);

confMat=confusionmat(refVec,cellVec);
totalPix=sum(sum(confMat));
confMatper(:,:,indx)=confMat./totalPix;
save(['/Volumes/Briggs_10TB/CRISPdata/Analysis/' title 'Voronoi_results.mat']);

 
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


