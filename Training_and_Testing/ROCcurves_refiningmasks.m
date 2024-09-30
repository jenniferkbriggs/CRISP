%%This script is how we define the training data and develop the ROC curve
close all
%clear all
clc

addpath('~/Documents/GitHub/UniversalCode/');
addpath('~/Documents/GitHub/Islet_Heterogeneity/')
%chose islets to analyze
filename = ["three","sample", "five", "two","one"];
csvname = ["H2BmCherry Ucn3GCaMP-3_Detailed.csv",  "H2BmCherry Ucn3GCaMP sample_Detailed.csv", "H2BmCherry Ucn3GCaMP-5_Detailed.csv", "H2BmCherry Ucn3GCaMP-2_Detailed.csv", "H2BmCherry Ucn3GCaMP-1_Detailed.csv"]
datapath = ['/Volumes/Briggs_10TB/Merrin/Confocal/'] 
savepath = ['~/Documents/GitHub/ST_Analysis/Data/']



%initialize which time index to use nuclear locations
%timetouse = 243%three: 217 - somewhat wiggly. %sample: 194; %five: 243%two: 363; %one: 193
timetouse = [217, 194, 243, 363, 193];
%Number of cells per islet for the training case
perislet = 15;

%set seed
thr = [-2:.05:3];%[0.25:0.05:1]

for l = [1:length(thr)]
for kt = 1:length(filename)

%set random number generator
rng(kt*2)

%First load good mask: 
load(strrep(strjoin([savepath filename(kt) '_Good.mat']), ' ', ''))
close all

reference_hold = CellMask;

for mm = [1:3]; %I will also compare doing st analysis on good mask
    reference = reference_hold;
    ca_files = dir(strrep(strjoin([datapath filename(kt) '/' '*C1*.tif']),' ',''));
    nuc_files = dir(strrep(strjoin([datapath filename(kt) '/' '*C2*.tif']),' ',''));

   
    if contains(ca_files(1).name, '._')
        ca_files(1) = [];
    end
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
    
    load(strrep(strjoin([savepath filename(kt) '_' masktype '.mat']), ' ', ''))
    

    %Opts: 
    Opts.fig = 0;
    Opts.st_thr = thr(l);
    Opts.Thr = 'st'
    %ST analysis:
    CellMask_st = STanalysis_refinemasks(Islet_vid, CellMask, Opts);

    %these pixels were in the reference mask but not in the mask that is
    %going through st analysis - therefore they cannot be removed so we
    %removed them from the percentages
    picsinreference_notinmask = (find(CellMask - reference < 0));

    if length(picsinreference_notinmask) > 0
        CellMask(picsinreference_notinmask) = [];
        CellMask_st(picsinreference_notinmask) = [];
        reference(picsinreference_notinmask)=[];
    end

    picsNotinMaskorReference = intersect(intersect(find(CellMask == 0), find(reference == 0)), find(CellMask_st == 0));
    
    CellMask(picsNotinMaskorReference)=[];
    CellMask_st(picsNotinMaskorReference)=[];
    reference(picsNotinMaskorReference)=[];

    pixremoved = find((CellMask_st-CellMask)~=0);
    %pix incorrectly classified should be removed: 
    incorrect = find(CellMask ~= 0 & reference~= 0 & CellMask~=reference);
    
    pix_shouldberemoved = find((CellMask-reference) > 0);
    pix_shouldberemoved = [pix_shouldberemoved, incorrect];
    pixnotremoved = intersect(find(CellMask), find(CellMask_st));

    pixshouldnotberemoved = setdiff(find(reference~=0), incorrect);
    

    %negative = remove, positive = keep

    %false negative: -> number of pixels removed that shouldn't have been
    fn = length(setdiff(pixremoved, pix_shouldberemoved));

    %true negative: -> percent of pixels correctly removed
    tn = length(intersect(pixremoved, pix_shouldberemoved));

    %false positive: -> percent of pixels not removed that should have been
    fp = length(intersect(pix_shouldberemoved, pixnotremoved));

    %true positive: -> pixels not removed that shoulve not have
    %been
    tp = length(intersect(pixnotremoved, pixshouldnotberemoved));

    tn_all(mm,l,kt) = tn;
    tp_all(mm,l,kt) = tp;
    fn_all(mm,l,kt) = fn;
    fp_all(mm,l,kt) = fp;

    sensitivity(mm, l, kt) = tp/(tp+fn);
    specificity(mm, l, kt) = tn/(tn+fp);

        disp(tp/(tp+fn)+tn/(tn+fp))
     if tp+fp+fn+tn ~= length(reference)
         keyboard
     end

    end
end
end

precision = tp_all./(tp_all+fp_all);
recall = sensitivity;


%turns out i defined positive and negative backward: 


figure, plot(squeeze(recall(1,:,:)), squeeze(precision(1,:,:)))

figure,
for i = [1,2,3,4,5]
    hold on, plot((squeeze(1-specificity(1,1:end,i))), (squeeze((sensitivity(1,1:end,i)))))
    disp(trapz((squeeze(1-specificity(1,1:end,i))), (squeeze((sensitivity(1,1:end,i))))))
end
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Bad Masks')
 set(gcf, 'color','white')
set(gca, 'box','off')
set(gca, 'FontSize',15)




figure,
for i = [1,2,3,4,5]
    hold on, plot((squeeze(1-specificity(2,1:end,i))), (squeeze((sensitivity(2,1:end,i)))))
    disp(trapz((squeeze(1-specificity(2,1:end,i))), (squeeze((sensitivity(2,1:end,i))))))
end
figure, plot(squeeze(1-specificity(2,:,:)), squeeze(sensitivity(2,:,:)))
xlabel('False Positive Rate')
ylabel('True Positive Rate')
title('Medium Masks')
 set(gcf, 'color','white')
set(gca, 'box','off')
set(gca, 'FontSize',15)


%% plot means: 
% Calculate mean and standard deviation
mean_x = mean(squeeze(1-specificity(2,:,:))');
std_x = std(squeeze(1-specificity(2,:,:))');
mean_y = mean(squeeze(sensitivity(2,:,:))');
std_y = std(squeeze(sensitivity(2,:,:))');

% Create the shaded area for the mean and standard deviation
figure;
hold on;

l3 = plot(squeeze(1-specificity(2,:,:)), squeeze(sensitivity(2,:,:)), 'k')

% Plot means
l1 =plot(mean_x(1:end-1), mean_y(1:end-1), 'k', 'LineWidth', 5);

% Customize plot
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('Medium Masks');
set(gcf, 'color', 'white');
set(gca, 'box', 'off');
set(gca, 'FontSize', 15);
ylim([0,1])
xlim([0,1])
hold on;



%plot means: 
% Calculate mean and standard deviation
mean_x = mean(squeeze(1-specificity(1,:,:))');
std_x = std(squeeze(1-specificity(1,:,:))');
mean_y = mean(squeeze(sensitivity(1,:,:))');
std_y = std(squeeze(sensitivity(1,:,:))');
% Create the shaded area for the mean and standard deviation
hold on;

% Plot means
l4 = plot(squeeze(1-specificity(1,:,:)), squeeze(sensitivity(1,:,:)), 'r')

l2 = plot(mean_x(1:end-1), mean_y(1:end-1), 'r', 'LineWidth', 5);

% Customize plot
xlabel('False Positive Rate');
ylabel('True Positive Rate');
set(gcf, 'color', 'white');
set(gca, 'box', 'off');
set(gca, 'FontSize', 20);
ylim([0,1])
xlim([0,1])
legend([l1, l2, l3(1), l4(1)], 'Medium Masks Mean','Bad Masks Mean', 'Individual Islets', 'Individual Islets')
title('Refining Masks: Receiver Operating Characteristic Curve')
hold off;

%% AUC: 
for i = 1:3
    for j = 1:2
        x = squeeze(1 - specificity(j, :, i));
        y = squeeze(sensitivity(j, :, i));
        
        % Ensure x is sorted in ascending order
        [x, sortIdx] = sort(x);
        y = y(sortIdx);
        
        % Calculate AUC using trapz
        auc(i, j) = trapz(x, y);
        
        % Check if AUC is negative
        if auc(i, j) < 0
            warning('AUC is negative for i=%d, j=%d. Please check the data.', i, j);
        end
    end
end

disp(['Mean AUC for bad masks: ' num2str(mean(auc(:,1))), ' with standard deviation: ' num2str(std(auc(:,1)))])
disp(['Mean AUC for medium masks: ' num2str(mean(auc(:,2))), ' with standard deviation: ' num2str(std(auc(:,2)))])


%% 

save Data/CorrelatioThresholdROC_st.mat sensitivity specificity thr tn_all fp_all tp_all fn_all






%after this also do threshold for correlation rather than standard
%deviation



