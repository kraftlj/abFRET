%%% Automated analysis of acceptor photobleaching
%%
%     abFRET Copyright (C) 2014  Lewis J. Kraft
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or (at
%     your option) any later version.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Running the script launches the acceptor phtobleaching quantification
% script.  It will prompt the user for the location of directory harboring
% the imaging files.  The script assumes there is a donor and an acceptor
% channel.
%
%
%% This is a script to analyze acceptor photobleaching data
%Inputs are images of donor, and images of acceptor both before and after
%bleaching.
%Outputs are:
% FRET efficiency images;
% % acceptor photobleached;
% E averaged over the whole cell;
% E averaged over the spots in the donor channel;
% E averaged over the spots in the acceptor channel;
%% Begin
clear all; close all; clc;
tic %start the clock
location = uigetdir;
files = dir(fullfile(location,'*.lsm')); %This catalogs all of the .lsm files in the location directory
if exist(fullfile(location,'SaveFolder')) %If the SaveFolder for the analysis exists remove it
    rmdir(fullfile(location,'SaveFolder'),'s');
end
mkdir(fullfile(location,'SaveFolder')); %Make a new SaveFolder
%% Image importing
%This catalogs all of the .lsm files in the location directory
wb=waitbar(0,'Overall Time Left'); %This initializes a waitbar
savedata = cell(length(files),5);
for index = 1:length(files); %This goes through the image files one by one
    % savedata = cell(2,5);
    % for index = 1:2; %This goes through the image files one by one
    wb=waitbar(index./length(files),wb); %This updates the waitbar each time it moves to the next image file
    clearvars('-except','savedata','location','files','fid','index','wb')
    data=bfopen(fullfile(location,files(index).name));
    %         for j=1:length(files)
    %         [num2str(j),' ', files(j).name]
    %         end
    % data=bfopen(fullfile(location,files(65).name));
    CFPprebleach=data{1,1}{9};
    YFPprebleach=data{1,1}{10};
    CFPpostbleach=data{1,1}{11};
    YFPpostbleach=data{1,1}{12};
    %     YFPsub=data{1,1}{30};
    
    clearvars('data')
    
    C=xcorr2(double(CFPpostbleach),double(CFPprebleach)); % Calculate the 2d cross correlation
    [max_cc, imax] = max(abs(C(:)));  % Find the maximum correlation
    [ypeak, xpeak] = ind2sub(size(C),imax(1));  % convert to x,y coordinates
    corr_offset = [ (ypeak-size(CFPprebleach,1)) (xpeak-size(CFPprebleach,2)) ];  % calculate the template offset
    I = CFPprebleach; % Load template
    tx = corr_offset(2); % Integer translation along the horizontal component
    ty = corr_offset(1); % Integer translation along the vertical component
    % Transform image
    t = maketform('affine',[1 0 ; 0 1; tx ty]);
    CFPprebleach = imtransform(I,t,'XData',[1 size(I,2)],'YData',[1 size(I,1)]); % Offset the CFPprebleach to match CFPpostbleach
    clearvars('C','max_cc','imax','ypeak','xpeak','corr_offset','I','t')
    
    h=figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,2,1)
    imagesc(YFPpostbleach)
    title('Draw the ROI on the image to the right')
    axis square
    subplot(1,2,2)
    imagesc(YFPprebleach)
    axis square
    BWuser=roipoly(YFPprebleach./(max(double(YFPprebleach(:)))/2^16));
    close(h);
    
    
    thresholding = 1;
    while thresholding == 1
        
        prompt = {'Enter background threshold:','Enter CFP spot threshold:','Enter YFP spot threshold'};
        dlg_title = 'Input';
        num_lines = 1;
        def = {'150','500','500'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        try
            
            UEB = medfilt2(CFPprebleach, [15 15],'symmetric');
            CFPspots = CFPprebleach-UEB;
            BWCFPspots = im2bw(uint16(CFPspots), str2double(answer{2})/2^16);
            BWCFPspots = imopen(BWCFPspots,strel('disk',2));
            
            UEB = medfilt2(YFPprebleach, [15 15],'symmetric');
            YFPspots = YFPprebleach-UEB;
            BWYFPspots = im2bw(uint16(YFPspots), str2double(answer{3})/2^16);
            BWYFPspots = imopen(BWYFPspots,strel('disk',2));
            
            BW=im2bw(UEB,str2double(answer{1})/2^16);
            BWbackground = ~imfill( bwareaopen(BW,4000),'holes');
            BW=~BWbackground;
            BW(~BWuser)=0;
            
        catch err
            
            if strcmp(err.identifier, 'images:im2bw:outOfRangeThreshLuminanceLevel')
                
                prompt = {'Enter background threshold:','Enter CFP spot threshold:','Enter YFP spot threshold'};
                dlg_title = 'Input';
                num_lines = 1;
                def = {'50','300','300'};
                answer = inputdlg(prompt,dlg_title,num_lines,def);
            end
        end
        
        
        
        clearvars('UEB','CFPspots','YFPspots')
        
        BWCFPspots(~BW)=0;
        BWCFPspots = imdilate(BWCFPspots,strel('disk',2));
        BWYFPspots(~BW)=0;
        BWYFPspots = imdilate(BWYFPspots,strel('disk',2));
        BW(BWCFPspots)=0;
        BW(BWYFPspots)=0;
        
        BWRGB(:,:,1)=uint8(BWYFPspots).*255;
        BWRGB(:,:,2)=uint8(BW).*255;
        BWRGB(:,:,3)=uint8(BWCFPspots).*255;
        
        h = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(2,2,1)
        imagesc(CFPprebleach)
        colormap('gray')
        title('CFPprebleach')
        axis image
        
        subplot(2,2,2)
        imagesc(YFPprebleach)
        colormap('gray')
        title('YFPprebleach')
        axis image
        
        subplot(2,2,3)
        imshow(BWRGB)
        title('Masks')
        axis image
        
        subplot(2,2,4)
        axis image
        texth=text(0.2,0.5,{files(index).name});
        set(texth,'interpreter','none');
        
        
        
        
        
        
        
        % Construct a questdlg with three options
        choice = questdlg('Are you satisfied with the thresholds?','Continue?','Yes', 'No','Yes');
        
        % Handle response
        switch choice
            case 'Yes'
                
                thresholding = 0;
                close (h);
            case 'No'
                
                thresholding = 1;
                close (h);
        end
        
    end
    
    BGCFPpre = mean(CFPprebleach(BWbackground));
    BGCFPpost = mean(CFPpostbleach(BWbackground));
    BGYFPpre = mean(YFPprebleach(BWbackground));
    BGYFPpost = mean(YFPpostbleach(BWbackground));
    
    CFPprebleach=double(CFPprebleach)-BGCFPpre;
    YFPprebleach=double(YFPprebleach)-BGYFPpre;
    CFPpostbleach=double(CFPpostbleach)-BGCFPpost;
    YFPpostbleach=double(YFPpostbleach)-BGYFPpost;
    
    thCFPpre = CFPprebleach;
    thCFPpost = CFPpostbleach;
    thYFPpre = YFPprebleach;
    thYFPpost = YFPpostbleach;
    thCFPpre(BWbackground) = nan;
    thCFPpost(BWbackground) = nan;
    thYFPpre(BWbackground) = nan;
    thYFPpost(BWbackground) = nan;
    binthCFPpre = imresize(thCFPpre,1/5);
    binthCFPpost = imresize(thCFPpost,1/5);
    ratioimg = (double(binthCFPpost)-double(binthCFPpre))./double(binthCFPpost);
    
    ratioimg(ratioimg(:)<0)=nan;
    ratioimg(ratioimg(:)>1)=nan;
    
    fretimage=ratioimg;
    fretimage(~imresize(BWuser,1/5))=nan;
    
    imwrite(uint8(fretimage.*255),fullfile(location,'SaveFolder',[files(index).name(1:end-4),'_FRET.tif']))
    imwrite(BWRGB,fullfile(location,'SaveFolder',[files(index).name(1:end-4),'_Masks.tiff']))
    imwrite(uint8(CFPprebleach.*(255/2^12)),fullfile(location,'SaveFolder',[files(index).name(1:end-4),'_CFPpre.tiff']))
    imwrite(uint8(CFPpostbleach.*(255/2^12)),fullfile(location,'SaveFolder',[files(index).name(1:end-4),'_CFPpost.tiff']))
    imwrite(uint8(YFPprebleach.*(255/2^12)),fullfile(location,'SaveFolder',[files(index).name(1:end-4),'_YFPpre.tiff']))
    imwrite(uint8(YFPpostbleach.*(255/2^12)),fullfile(location,'SaveFolder',[files(index).name(1:end-4),'_YFPpost.tiff']))
    
    Accpost = mean(YFPpostbleach(BWuser))-BGYFPpost;
    Accpre = mean(YFPprebleach(BWuser))-BGYFPpre;
    Donpost = mean(CFPpostbleach(BW))-BGCFPpost;
    Donpre = mean(CFPprebleach(BW))-BGCFPpre;
    DonpostCFPspot = mean(CFPpostbleach(BWCFPspots))-BGCFPpost;
    DonpreCFPspot = mean(CFPprebleach(BWCFPspots))-BGCFPpre;
    DonpostYFPspot = mean(CFPpostbleach(BWYFPspots))-BGCFPpost;
    DonpreYFPspot = mean(CFPprebleach(BWYFPspots))-BGCFPpre;
    %Calculation for Percent Acceptor Bleaching
    PAB = 100*(Accpre - Accpost)./Accpre;
    %Calculation for Energy Transfer
    ETcyt = 100*(Donpost - Donpre)./Donpost;
    ETCFPspots = 100*(DonpostCFPspot-DonpreCFPspot)./DonpostCFPspot;
    ETYFPspots = 100*(DonpostYFPspot-DonpreYFPspot)./DonpostYFPspot;
    
    savedata(index,:)={files(index).name,PAB,ETcyt,ETCFPspots,ETYFPspots};
    
    
end

savedata=savedata';
header={'FileNames','Percent acceptor bleached','Cytosol E','CFP spots E','YFP spots E'};
fid = fopen(fullfile(location,'Results.txt'),'w');
fprintf(fid, '%s\t %s\t %s\t %s\t %s\r\n', header{:});
fprintf(fid, '%s\t %g\t %g\t %g\t %g\r\n', savedata{:});
fclose(fid);

close(wb);

toc
