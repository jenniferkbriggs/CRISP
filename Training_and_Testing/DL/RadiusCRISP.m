
close all
clearvars -except radiusoutput
clc

%Find the radial location of training
th_pix = 0.8183; %threshold for pixels
th_rad = 0.83; %score for radius. 

%load:
savepath = '/Users/brigjenn/Documents/GitHub/ST_Analysis/Data/'
load(([savepath 'sample_Good.mat']))

kl = 4
ll =14
files = dir('/Volumes/Briggs_10TB/CRISPdata/Analysis/*.csv')
CellPose = load(['/Volumes/Briggs_10TB/CRISPdata/Analysis/' files(ll).name 'MaskWithReferenceCells.mat']);


% need to match reference cells with real cell - loop over to find the
% cells with the largest overlap first:
    CellPoseMask_updated = zeros(size(CellMask));
usedindx = [];
for i = 1:max(unique(CellMask))
    ref = zeros(size(CellMask));
    [x,y]=find(CellMask==i);
    for j = 1:length(x)
    ref(x(j),y(j)) = 100;
    end


    diff = ref-(CellPose.CellMask); 
    overlap = diff(find(diff>5&diff<100));
    CellMaskIndx = 100-mode(overlap);


    if length(intersect(usedindx, CellMaskIndx))==0
    [x,y]=find(CellPose.CellMask==CellMaskIndx);
    for j = 1:length(x)
    CellPoseMask_updated(x(j),y(j)) = i;
    end
    end
    usedindx=[usedindx, CellMaskIndx]



    %look for positive numbers less than 100 and then find the one with the
    %highest count
end
% 
% 
% figure , nexttile, imshow(CellMask),    for c = 1:max(max(CellMask))
% [x,y] = find(CellMask ==c)
% text(mean(y),mean(x), num2str(c)); % Labels cells in the image with their respective region number
% end,
% nexttile, imshow(CellPoseMask_updated),   for c = 1:max(max(CellPoseMask_updated))
% [x,y] = find(CellPoseMask_updated ==c)
% text(mean(y),mean(x), num2str(c)); % Labels cells in the image with their respective region number
% end


%Two anonymous functions that calculate the pixels within the radius.
%Circfilled gives all pixels inside circle, Circ gives only pixels on the
%circumference
circfilled = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

ct = 1;



for j = 1:max(max(CellMask)) %loop over cells
    %cut islet video:

    [TrueCellx TrueCelly] = find(CellMask == j);        
    %find maximum radius:
    maxedout = 0;
    radius = 1;
    while maxedout == 0
        pixy = circfilled(radius, (fliplr(NucLoc_keep(j,:))));
        if length(intersect(pixy, [TrueCellx TrueCelly], 'rows')) ~= length(pixy) %then we've hit the maximum radius
            maxedout = 1;
        else
            radius = radius +1;
        end
    end
        trueRadii(j)= radius;

 %run CRISP
    ct = ct+1;
    opts.ca = ca_im;
    opts.gif = 0;
    opts.normalized = 1;
    opts.th_pix = th_pix;
    opts.score_thr = th_rad;
    opts.radiusstart = 2; %how large to make baseline radius
    opts.title = [filename(kt) 'Islet' num2str(i)]
    [Correlation, radius, Pixelsx, Pixelsy, score] = CRISP_annulus_corr((Islet_vid),(NucLoc_keep(j,:)),opts);
    Corr_all(j).corr = Correlation;
    %CRISP_rad(j)=radius;
    CRISP_rad_all(j,:) = [3:42];
    CRISP_score_all(j,:) = score;

    %radius of DL:
    [TrueCellx TrueCelly] = find(CellPoseMask_updated == j); 
    if isempty(TrueCellx)
        DL_rad(j) = 0;
    else
    %find maximum radius:
    maxedout = 0;
    radius = 1;
    while maxedout == 0
        pixy = circfilled(radius, (fliplr(NucLoc_keep(j,:))));
        if length(intersect(pixy, [TrueCellx TrueCelly], 'rows')) ~= length(pixy) %then we've hit the maximum radius
            maxedout = 1;
        else
            radius = radius +1;
        end
    end
        DL_rad(j)= radius;
    end
end

ScoresinRadius = nan(size(CRISP_score_all));
ScoresoutofRadius = nan(size(CRISP_score_all));
for j = 1:size(CRISP_score_all,1)
    ScoresinRadius(j,1:trueRadii(j)-2) = CRISP_score_all(j,1:trueRadii(j)-2);
    ScoresoutofRadius(j,trueRadii(j)-1:end) = CRISP_score_all(j,trueRadii(j)-1:end);
end


%true positive is radius is smaller than true radius:
ct=1
for i = [0:0.01:1]
    tp(ct) = length(find(ScoresinRadius > i))./length(find(~isnan(ScoresinRadius)));
    fp(ct) = length(find(ScoresoutofRadius > i))./length(find(~isnan(ScoresoutofRadius)));
    tn(ct) = length(find(ScoresoutofRadius < i))./length(find(~isnan(ScoresinRadius)));
    fn(ct) = length(find(ScoresinRadius < i))./length(find(~isnan(ScoresoutofRadius)));
    ct = ct+1;
end

radiusoutput.fn(kl,:) = fn;
radiusoutput.tp(kl,:) = tp;
radiusoutput.tn(kl,:) = tn;
radiusoutput.fp(kl,:) = fp;
radiusoutput.referencescores(kl,:) =  [0:0.01:1];
radiusoutput.DL_rad_all(kl, :) = DL_rad;
radiusoutput.True_rad_all(kl, :) = trueRadii;
%output.CRISP_rad_all(kl, :) = CRISP_rad;
radiusoutput.missingDLcells(kl) = length(find(DL_rad == 0))
radiusoutput.CRISP_score_all(:,:,kl) = CRISP_score_all;

save(['/Volumes/Briggs_10TB/CRISPdata/Analysis/' files(ll).name 'RadiusAnalysis.mat']);

%% 
th_rad = 0.83
for i = 1:5
    for j = 1:15
        foo = find(radiusoutput.CRISP_score_all(j,:,i) < th_rad);
        CRISP_rad(j,i) = foo(1)+3
    end
end

mean2(CRISP_rad - radiusoutput.True_rad_all')
mean2(radiusoutput.DL_rad_all' - radiusoutput.True_rad_all')

std(CRISP_rad - radiusoutput.True_rad_all', [], 'all')
std(radiusoutput.DL_rad_all' - radiusoutput.True_rad_all', [], 'all')
