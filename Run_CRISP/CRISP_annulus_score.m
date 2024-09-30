function [radii] = CRISP_annulus_corr(images,NucLoc,opts)
%Jennifer Briggs 03.2022
% Here the score (p+) is calculated and the semi-minor axis length is identified

radii = 1;
score = 1;
%set threshold
while score > th_rad
    %if percent of pixels in radius is greater than threshold, we can increase radius
     [pix_belowth(radii)] = length(find(nonzeros(Correlation(:,radii)) > th_pix));
     [tot_pix] = length(nonzeros(Correlation(:,radii)));
     score = pix_belowth(radii)/tot_pix-radii/100;
     radii = radii+1;
end

end