function [] = gen_Baseline_For_DFOF(  )
%%gen_Baseline_For_DFOF generate the base line image by averaging all the time points. 
% This is needed for delta F over F calculations.
% Author: Raju Tomer (rajutomer@gmail.com)
%%

%% Parameters %%
dir_base = 'D:\SPED_data\Exp1'; % base dir
dir_data = [dir_base '\LOC000_dc']; % input data dir
dir_out = [dir_base '\ref_im']; % output dir
stack_size = [2048 632 40]; % stack size
DC_IT = 10; % used for parsing files which were generated by using particular iterations of RL deconvolutions.
if_tiff_stack = 1; % if files are in 3d tiff stack format
if_reverse_stack = 0; % if to reverse the slice order in the stack
if_raw_data = 0; % if processing the raw data or the deconvolved data, needed for file name parsing
if_med_filt = 0; % if to use median filtering
start_image_ind = 1; % start image ind (coressponding to time point)
end_image_ind = 4400; % end image ind
%%

if (if_raw_data == 1)
    list = dir([dir_data '\*.stack']);
else
    list = dir([dir_data '\DC' num2str(DC_IT) '*.tif']);
end

avg = [];
init_sw = 0;
num_images = end_image_ind - start_image_ind + 1;

for i = start_image_ind:end_image_ind
    fname = [dir_data '\' list(i).name];
    im = [];
    if (if_tiff_stack > 0)
        im = zeros([stack_size(2) stack_size(1) stack_size(3)], 'uint16');
        for r = 1:stack_size(3)
            im(:,:,r) = imread(fname,r);
        end
    else
        fid = fopen(fname, 'r', 'l' );
        im = fread( fid, prod(stack_size), 'uint16' );
        fclose( fid );
        im = reshape(im,stack_size);
        im = permute(im, [2 1 3]);    
    end
    if (if_reverse_stack == 1)
        im = flip(im,3);
    end
    if (if_med_filt == 1)
        for rj = 1:stack_size(3)
            im(:,:,rj) = medfilt2(im(:,:,rj));
        end
    end

    if (init_sw == 0)
        avg = double(im)./num_images;
        init_sw=1;
    else
        avg = avg + double(im)./num_images;
    end
end

ofname = '';
if (if_raw_data == 1)
    if (if_med_filt == 1)
        ofname = [dir_out '\Avg_Raw_medF' '_Range' num2str(start_image_ind) '-' num2str(end_image_ind) '.tif'];
    else
        ofname = [dir_out '\Avg_Raw' '_Range' num2str(start_image_ind) '-' num2str(end_image_ind) '.tif'];
    end
else
    if (if_med_filt == 1)
        ofname = [dir_out '\Avg_medF_DC' num2str(DC_IT) '_Range' num2str(start_image_ind) '-' num2str(end_image_ind) '.tif'];
    else
        ofname = [dir_out '\Avg_DC' num2str(DC_IT) '_Range' num2str(start_image_ind) '-' num2str(end_image_ind) '.tif'];
    end
end

for k = 1:stack_size(3)
 	imwrite(uint16(avg(:,:,k)), ofname, 'tiff', 'Compression', 'none', 'WriteMode', 'append');
end