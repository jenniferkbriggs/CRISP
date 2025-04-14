
close all
clearvars -except radiusoutput
clc

%Find the radial location of training
th_pix = 0.8183; %threshold for pixels
th_rad = 0.83; %score for radius. 

%load:
savepath = '/Users/brigjenn/Documents/GitHub/ST_Analysis/Data/'
load(([savepath 'five_Good.mat']))

ll = 5

%files = dir('/Volumes/Briggs_10TB/CRISPdata/Analysis/*.csv')

%load voronoi:%
%Voron = load(['/Volumes/Briggs_10TB/CRISPdata/Analysis/' title 'Voronoi_results.mat']);
files = dir('/Users/brigjenn/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Anschutz/Islet/CRISP/ReviewFigs/SelectedNucs/*.csv')
cMask = load([files(1).folder '/' files(ll).name])
Voron.CellMask = cMask;


% need to match reference cells with real cell - loop over to find the
% cells with the largest overlap first:
    VoronMask_updated = zeros(size(CellMask));
usedindx = [];
for i = 1:max(unique(CellMask))
    ref = zeros(size(CellMask));
    [x,y]=find(CellMask==i);
    for j = 1:length(x)
    ref(x(j),y(j)) = 100;
    end


    diff = ref-(Voron.CellMask); 
    overlap = diff(find(diff>5&diff<100));
    CellMaskIndx = 100-mode(overlap);


    if length(intersect(usedindx, CellMaskIndx))==0
    [x,y]=find(Voron.CellMask==CellMaskIndx);
    for j = 1:length(x)
    VoronMask_updated(x(j),y(j)) = i;
    end
    end
    usedindx=[usedindx, CellMaskIndx];



    %look for positive numbers less than 100 and then find the one with the
    %highest count
end
% 
% 
% figure , nexttile, imshow(CellMask),    for c = 1:max(max(CellMask))
% [x,y] = find(CellMask ==c)
% text(mean(y),mean(x), num2str(c)); % Labels cells in the image with their respective region number
% end,
% nexttile, imshow(VoronMask_updated),   for c = 1:max(max(VoronMask_updated))
% [x,y] = find(VoronMask_updated ==c)
% text(mean(y),mean(x), num2str(c)); % Labels cells in the image with their respective region number
% end


%Two anonymous functions that calculate the pixels within the radius.
%Circfilled gives all pixels inside circle, Circ gives only pixels on the
%circumference
circfilled = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');
circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)],'rows');

ct = 1;



for j = 1:max(max(CellMask)) %loop over cells

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

    %radius of DL:
    [TrueCellx TrueCelly] = find(VoronMask_updated == j); 
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

mean2(DL_rad - trueRadii')

std(DL_rad - trueRadii, [], 'all')


save(['/Volumes/Briggs_10TB/CRISPdata/Analysis/' files(ll).name 'RadiusAnalysis.mat']);
