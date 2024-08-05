function [CellMask_updated] = STDanalysis_refinemaks(images, CellMask, Opts)
%Jennifer Briggs 02.2022
%This function takes a time course of pixels with masks differentiating
%cells and removes pixels within that mask that do not correlate with the
%rest of the pixels (e.g. likely not a part of the cell.)

%images is the image file in matrix form (x,y,t). Currently only accepting
%3D matrices (e.g. no z stacks) so z stacks should be fed in separately. 

%Cell Mask is a pixel x pixel matrix of the imaging file. Any pixel that is
%considered to be a part of a cell has the index of that corresponding
%cell. For example, if we are looking at a 5x5 pixel images with 1 cell
%inside, cell mask would be: 
% [0, 0, 0, 0, 0; 0, 1, 0, 0, 0; 0, 1, 1, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0,
% 0]


%Opts - Opt.fig = 0 or 1 for no figures or figures respectively
   %    Opt.st_thr = double that defines how many standard deviations to
   %    count as bad pixels. 
    images = double(images)+0.01; %note that it rotates 90 degrees again. not sure why
    CellMasksave = CellMask; %save old cell mask

    if Opts.fig
    [x,y] = find(CellMask)
    figure,plot(x,y,'o')
    end

     %gif('FullGif.gif')
    numcells = unique(CellMask);
    for i = 1:length(numcells)-1
        
        cctt = 1;
        badpix = [2 2]; %arbitrary array that will be filled with the bad pixel locations
        %Get video of the very center of the mask:

        TCMask = CellMask; %Pulls in CellMask array
        TCMask(TCMask ~= i) = 0; %Gets rid of all masks besides current one
        MaskedIMGstack = images.*logical(TCMask); %Applies current mask to all frames of timecourse
        
        for ii = 1:size(MaskedIMGstack,3)
            TCnoZero = MaskedIMGstack(:,:, ii); %Pulls in current frame
            TCnoZero = TCnoZero(TCnoZero>0); %Accounts for any zeros from preallocation
            if size(TCnoZero,1) > size(TCnoZero,2)
                TCcheck(ii,:) = TCnoZero';
            else
                TCcheck(ii,:) = TCnoZero;
            end
            %TCmiddle(ii) = mean2(Middleimage(:,:,ii));
            TC(ii) = mean(TCnoZero);
        end

        
       [rr, cc] = find(mean(MaskedIMGstack,3));

        xx = round(mean(rr)-std(rr)./2:round(mean(rr))+std(rr)./2);
        yy = round(mean(cc)-std(cc)./2:round(mean(cc))+std(cc)./2);

        %make grid:
        [xxx, yyy] = meshgrid(xx,yy);

        Middleimage = MaskedIMGstack(reshape(xxx,1,[]), ...
             reshape(yyy,1,[]),:);


        TCreference = (squeeze(mean(mean(Middleimage))));

        while ~isempty(badpix)
        TCMask = CellMask; %Pulls in CellMask array
        TCMask(TCMask ~= i) = 0; %Gets rid of all masks besides current one
        MaskedIMGstack = images.*logical(TCMask); %Applies current mask to all frames of timecourse
        
        for ii = 1:size(MaskedIMGstack,3)
            TCnoZero = MaskedIMGstack(:,:, ii); %Pulls in current frame
            TCnoZero = TCnoZero(TCnoZero>0); %Accounts for any zeros from preallocation
            if size(TCnoZero,1) > size(TCnoZero,2)
                TCcheck(ii,:) = TCnoZero';
            else
                TCcheck(ii,:) = TCnoZero;
            end
            %TCmiddle(ii) = mean2(Middleimage(:,:,ii));
            TC(ii) = mean(TCnoZero);
        end

       
        

        
        if size(TCcheck,2) >0
            try
                for c=1:size(TCcheck,2)
                        [C1(c)] = corr(TCcheck(:,c),TCreference);
                end
            catch
                keyboard
            end
        %pixel is bad if it shows no oscillations
        
        %badpix = find(Maxcor<5e6);
        if cctt == 1
            C1_references = Opts.st_thr;
        end
        %if length(badpix) > .75*size(TCcheck,2)    

                %standard deviation threshold - C1_references = std
               %badpix = find(C1<mean(C1) - C1_references);
             %Correlation threshold
               badpix = find(C1<C1_references);

               cctt = cctt +1;
               clear C1
        
        
        
       %end
       if cctt < 25

        %title(['Itteration' num2str(cctt)]), plot(rr,cc,'o')
        try
        for bb= 1:length(badpix)
            if Opts.fig
            hold on, plot(rr(badpix(bb)),cc(badpix(bb)),'ro','markerfacecolor','r')
            %saveas(gcf, ['Itteration' num2str(cctt) '.png']) 
            drawnow
            %gif
            end
            CellMask(rr(badpix(bb)),cc(badpix(bb))) = 0;
        end
        catch
            keyboard
        end
     
       else
    %        disp('Bad Pixels are too many')
    %        disp(savepath)
    %        fig = figure, plot(TCchecksm)
    %        x = input('Keep deleting or move on? 1 for move on, 0 for keep deleting')
    %        close(fig)
    %         if x == 0
    %             cctt = 1;
    %         else
    %             cctt = 5;
                badpix = [];
    %         end
       end
            clear TCcheck Maxcor Maxcor2

            %figure, plot(TCchecksm)
         else
           disp('No pixels')
           badpix = [];
           clear TCcheck
        end
        end


        CellTC(:,i) = TC; %Updates CellTC array with mean intensities from each frame of each cell
        %PlotLabels{i} = ['Cell' num2str(i)]; %Updates labels with current cell number for legend
    end
    CellMask_updated = CellMask;
end
% 
% correlation_data = C1;
% % Normalize the correlation data to range [0, 1]
% normalized_data = (correlation_data - min(correlation_data)) / (max(correlation_data) - min(correlation_data));
% 
% % Define a colormap (e.g., 'jet')
% colormap_name = 'jet';
% colormap_values = colormap(colormap_name);
% 
% % Map the normalized data to colormap indices
% colormap_indices = round(normalized_data * (colormap_length - 1)) + 1;
% 
% % Get the RGB values from the colormap
% rgb_values = colormap_values(colormap_indices, :);
% 
% % Display the RGB values
% disp(rgb_values);
