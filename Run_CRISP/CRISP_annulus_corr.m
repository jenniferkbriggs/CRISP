function [Correlation, radius, pixelsx, pixelsy] = CRISP_annulus_corr(images, NucLoc, opts)
    % Jennifer Briggs 03.2022
    % This function identifies the cell boundary based on pixel behavior. It expands a radius around the nucleus location and calculates correlations between pixel intensity changes, stopping when the behavior deviates beyond a threshold.
    
    % Parameters:
    % images: 3D matrix (x, y, t) of image data over time (no z-stacks).
    % NucLoc: [X,Y] position of the nucleus center. Median value is suggested if it changes over time.
    % opts: structure with additional options.
    % - opts.figs: 1 to generate a gif visualizing pixel removal.
    % - opts.score_thr: threshold for stopping the radius expansion.
    % - opts.th_pix: pixel correlation threshold for stopping.
    
        % Threshold for pixel correlation and score
        th_pix = opts.th_pix;   
        score_thr = opts.score_thr;
    
        % Initialize radius and maximum allowed radius
        radius = 2; % Initial radius
        maxradius = 40; % Maximum radius to avoid excessive expansion
    
        g = 0; % GIF counter for optional GIF creation
    
        % Pull in baseline oscillation data at the initial radius
        [image_base, allpix] = getcal_radius(radius, images, NucLoc, 1);
    
        try % Optionally create a GIF of the expanding nucleus region
            if opts.gif
                disp([opts.title '.gif'])
                giffig = figure;
                ti = tiledlayout(2,2); ti.TileSpacing = 'compact';
                
                % Mark and visualize pixels in the initial radius
                caim = insertMarker(opts.ca, allpix, 'color', 'r', 'Size', 1);
                nexttile(1), imshow(caim)
                xlim([NucLoc(1)-maxradius, NucLoc(1)+maxradius])
                ylim([NucLoc(2)-maxradius, NucLoc(2)+maxradius])
                
                % Plot the baseline oscillation in the next tile
                nexttile(2), plot([1:size(image_base, 1)], image_base)
                sgtitle([opts.title])
                
                % Save the GIF
                gif(['/Volumes/Briggs_10TB/Merrin/Confocal/AllGif' opts.title '.gif'])
                g = 1;
            end
        end
        
        % Average pixel behavior over time
        Y = mean(image_base, 2);
        
        % Initialize variables for storing correlation and pixel positions
        p = 0; j = 1;
        check = 0;
        Correlation = zeros(maxradius^2, maxradius);
        pixelsx = zeros(maxradius^2, maxradius);
        pixelsy = zeros(maxradius^2, maxradius);
    
        % Start CRISP (Cell Radius Identification by Simulated Pixels):
        score_thr = 1;
        Correlation = [];
        
        while score < score_thr % Expand until error exceeds the threshold
            clear Correlation
            radius = radius + 1; % Increase the radius at each iteration
    
            % Get all pixels within the current radius
            [image_new, allpix_new] = getcal_radius(radius, images, NucLoc, 0);
            
            % Calculate the average pixel behavior within the new radius
            Y = mean(image_new, 2);
    
            % Compute the correlation of new pixels with the average
            [r, p] = corr(image_new, Y);
            Correlation(1:length(r)) = r;
            pixelsx(1:length(allpix_new), j) = allpix_new(:, 1); % Store X coordinates of pixels
            pixelsy(1:length(allpix_new), j) = allpix_new(:, 2); % Store Y coordinates of pixels
    
            % Calculate the score: p+
            pix_belowth = length(find(nonzeros(Correlation) > th_pix));
            [tot_pix] = length(nonzeros(Correlation));
            score = pix_belowth(radius) / tot_pix - radius / 100;
        end
    
    end
    
    % Helper function to calculate pixel fluoresence for a given radius
    function [casmooth, allpix] = getcal_radius(radius, images, NucLoc, filled)
        % This function finds pixels around the nucleus in a circle with a given radius and extracts calcium oscillations for those pixels and smooths the data.
    
        if filled % If filled, include all pixels within the radius, not just the boundary
            circ = @(radius, NucLoc) unique([reshape((round([1:radius]'.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round([1:radius]'.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)], 'rows');
        else
            circ = @(radius, NucLoc) unique([reshape((round(radius.*cos(0:pi/2000:2*pi)+NucLoc(1))),[],1), reshape(round(radius.*sin(0:pi/2000:2*pi)+NucLoc(2)),[],1)], 'rows');
        end
        
        % Get all pixels within the radius
        allpix = circ(radius, NucLoc);
    
        % Extract calcium oscillations for each pixel
        for i = 1:length(allpix)
            calcium(:, i) = squeeze(images(allpix(i, 1), allpix(i, 2), :));
        end
         
        % Smooth the calcium data over time using a Gaussian filter
        casmooth = smoothdata(calcium, 'gaussian'); 
    end
    