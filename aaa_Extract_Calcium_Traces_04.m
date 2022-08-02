%% Clear
clear; 

%% Initial Setup parameters

% This GUI replaces the manual parameter setting in the ".m" file 
aaa_Initialize_Parameters_02

% Clean the workspace before saving it
% clr_ws = true;            

% Save the workspace as a matfile
% save_mat = true;       % If "True" at the end the script will clean up the workspace and save the rest as a matfile
                       % with the same filename at the same directory

% Manualy refine the ROI selection 
% manual_refine = false; % If set "True" there will be a window opened to manally refine the 
                       % found cells for the next process. There will
                       % instruction to show how to add more cells and
                       % remove onwanted

% Crop the raw image to save time or avoid bad ROI detection
% crop_image = false;  % Crop the loaded image to avoid unwanted regions 
                     % or artifacts. The number of line after cropping 
                     % should be even for the bidirectional offset 
                     % correction to work.



% Motion correction 
% do_mtn_crcn = true; % do_mtn_crcn : Do Motion Correction

% Scan phase correction
% do_scn_phs_crcn = true; % do_scn_phs_crcn : do scan phase correction 

% Save a movie or not of the reconstructed image
% mk_mv = false;  

% Number of frames to average 
% n_avrg = 6;            

% Get the full directory to the ScanImage tif file 
[tp_filename, tp_path] = uigetfile('.tif', ['Select the tiff image ' ...
    'from ScanImage']);
tp_filename = strcat(tp_path, tp_filename);

% Do you want to break the segment?
% brk_stck = false; % brk_stck: Break the stack
%-- Set the number of segments you wish to break the recording to
% n_sgmt = 5; % Number of segments
%-- Set the desired segment to be analyzed
% trgt_sgmt = 1; % target segment

% Gaussian smoothing in three dimension. 
% gssn_xyz = [1 0.9 3];

% Image histogram stretching threshhold 
% When the pixel value count falls below this threshhold, the 
% value will be set as maximum value.
% hstm_stch_tshd = 2000; 

% Image histogram streching threshhold after guassian smoothing 
% gssn_hstm_stch_tshd = 5000;

gcp;                            % start cluster
clear ans

%% Reading the file, averaging and stretching the histogram

% Read the stack
if brk_stck
    frm_rt = 30; % Frame Rate
    rc_tm = 30; % Recording time
    strt_frm = (trgt_sgmt -1) * (rc_tm * 60 * frm_rt/n_sgmt) + 1;  % Starting Frame 
    % (recording time in minute * seconds in a minute * frame rate = total number of frames)
    n_frm = (rc_tm * 60 * frm_rt/n_sgmt); % Number of frames to read
    im_stck = read_file(tp_filename, strt_frm, n_frm);
else
    im_stck = read_file(tp_filename);
end

% Crop if required 
if crop_image 
    % Set the crop region                     
    x1 = crop_coordinate(1);
    x2 = crop_coordinate(2);
    y1 = crop_coordinate(3);  % one of the y1 or y2 should be odd and the other even
    y2 = crop_coordinate(4);  % so the final croped image has an even number of lines.
                              % Even lines are necessary for the bidirectional correction
    if mod(abs(y1 - y2), 2) == 0 
        y1 = y1 + 1; 
    end
    im_stck = im_stck(y1:y2, x1:x2, :); % Crop image to remove unwanted regions
end 

% Clean up negative values
im_stck = single(im_stck);
im_stck = im_stck - min(im_stck, [], 'all');
im_stck = im_stck/max(im_stck, [], 'all');

%averaging along time 
[d1,d2,T] = size(im_stck);
T = T/n_avrg;
im_stck = mean(reshape(im_stck, d1, d2, n_avrg, T), 3);
im_stck = reshape(im_stck, d1, d2, T);
im_stck = im_stck/n_avrg;
[d1,d2,T] = size(im_stck);                          % dimensions of dataset
d = d1*d2;                                          % total number of pixels

%Streching the histogram 
[hist_cnt, hist_lctn] = imhist(im_stck);
hist_cnt_low_ix = find(hist_cnt<hstm_stch_tshd);
[~, hist_mx_ix] = max(hist_cnt);
low_in = hist_lctn(max(hist_cnt_low_ix(hist_cnt_low_ix<hist_mx_ix)));
high_in = hist_lctn(min(hist_cnt_low_ix(hist_cnt_low_ix>hist_mx_ix)));

im_stck = imadjustn(im_stck, [low_in high_in], []);

clear high_in hist_cnt hist_cnt_low_ix hist_lctn hist_mx_ix low_in ...
    x1 x2 y1 y2
%% Corret the every other line offset due to bi-directional scanning
if do_scn_phs_crcn
    %     nFrames = 200;
    nFrames = phs_crcn_frm; % Phase correction number of frames
    bidir_us = 3; % up-sampling to reach a sub-pixel image registration
    tic;[col_shift, M] = correct_bidirectional_offset(im_stck,nFrames,bidir_us);toc
    fprintf('Offset %1.1d pixels due to bidirectional scanning detected. \n',col_shift);

    figure;
    subplot(1,2,1)
    imshow(mean(im_stck(:,:,1:nFrames),3),[])
    title('raw')
    subplot(1,2,2)
    imshow(mean(M(:,:,1:nFrames),3),[])
    title('phase corrected')

    im_stck = M;
    clearvars M bidir_us col_shift nFrames
end

%% Smoothing the data

im_stck = imgaussfilt3(im_stck,gssn_xyz);

[hist_cnt, hist_lctn] = imhist(im_stck);
hist_cnt_low_ix = find(hist_cnt<gssn_hstm_stch_tshd);
[~, hist_mx_ix] = max(hist_cnt);
low_in = hist_lctn(max(hist_cnt_low_ix(hist_cnt_low_ix<hist_mx_ix)));
if isempty(low_in)
    low_in = 0;
end
high_in = hist_lctn(min(hist_cnt_low_ix(hist_cnt_low_ix>hist_mx_ix)));

im_stck = imadjustn(im_stck, [low_in high_in], []);

im_stck = im2uint8(im_stck);

clear high_in hist_cnt hist_cnt_low_ix hist_lctn hist_mx_ix low_in
%% set parameters for motion correction (first try out rigid motion correction)
if do_mtn_crcn
    options_rigid = NoRMCorreSetParms('d1',size(im_stck,1),'d2',size(im_stck,2),...
        'bin_width',50,...
        'max_shift', [50, 50, 5],...
        'us_fac',2,...
        'init_batch',200,...
        'boundary', 'copy', ...
        'upd_template', 0);
end 
%% perform motion correction
if do_mtn_crcn
    tic;
    [M,shifts,~,options_rigid] = normcorre(im_stck,options_rigid);
    toc

    % Run motion correction for a second round
    % YY = M;
    % clear M
    % [M,~,~,options_rigid] = normcorre(YY,options_rigid);

    figure;
    subplot(1,2,1)
    imshow(max(im_stck,[],3),[])
    title('max(stack), before motion correct')
    subplot(1,2,2)
    imshow(max(M,[],3),[])
    title('max(stack), after motion correct')

    im_stck = M;

    clear M
end
%% Set parameters for cell unit detection

% K = 100;                                         % number of components to be found
K = n_nrn;
% tau =  5;                                       % std of gaussian kernel (half size of neuron) 
tau = nrn_half;
% p = 2;
p = ca_prfl;
nb = bg_cmpt; 

options = CNMFSetParms(...   
    'd1',d1,'d2',d2,...                         % dimensionality of the FOV
    'p',p,...                                   % order of AR dynamics  1: only decay 2:rise and decay  
    'gSig',tau,...                              % half size of neuron
    'merge_thr',0.80,...                        % merging threshold  
    'nb',nb,...                                  % number of background components    
    'min_SNR',3,...                             % minimum SNR threshold
    'space_thresh',0.5,...                      % space correlation threshold
    'cnn_thr',0.2,...                           % threshold for CNN classifier    
    'ssub', 1, ...                              %spatial downsampling factor (default 1). 
    'tsub', 1 ...                               %temporal downsampling factor (default 1). 
    );
%% Data pre-processing

[P,im_stck] = preprocess_data(im_stck,p); 

%% fast initialization of spatial components using greedyROI and HALS

im_stck = single(im_stck);

tic; 
[Ain,Cin,bin,fin,center] = initialize_components(im_stck,K,tau,options,P);
toc;  

% display centers of found components
Cn = max(im_stck,[],3);
figure;
imagesc(Cn);
axis equal; 
axis tight; 
hold all;
scatter(center(:,2),center(:,1),'mo');
title('Center of ROIs found from initialization algorithm');

%% manually refine components (optional)
refine_components = manual_refine;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(im_stck,Ain,Cin,center,Cn,tau,options);
end
    
%% update spatial components
Yr = reshape(im_stck,d,T);
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

%% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% classify components

tic; 
rval_space = classify_comp_corr(im_stck,A,C,b,f,options); 
toc;

ind_corr = rval_space > options.space_thresh;   % components that pass the correlation test
                                                % this test will keep processes
                                        
%% further classification with cnn_classifier
try  % matlab 2017b or later is needed
    [ind_cnn,value] = cnn_classifier(A,[d1,d2],'cnn_model',options.cnn_thr);
catch
    ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
end     
                            
%% event exceptionality

fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
ind_exc = (fitness < options.min_fitness);

%% select components

keep = (ind_corr | ind_cnn) & ind_exc;

%% display kept and discarded components
A_keep = A(:,keep);
C_keep = C(keep,:);
figure;
subplot(121); 
montage(extract_patch(A(:,keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15]);
title('Kept Components');
subplot(122); 
montage(extract_patch(A(:,~keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15])
title('Discarded Components');
%% merge found components
[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A_keep,b,C_keep,f,P,S,options);

%%
display_merging = 1; % flag for displaying merging example
if and(display_merging, ~isempty(merged_ROIs))
    i = 1; %randi(length(merged_ROIs));
    ln = length(merged_ROIs{i});
    figure;
    set(gcf,'Position',[300,300,(ln+2)*300,300]);
    for j = 1:ln
        subplot(1,ln+2,j); imagesc(reshape(A_keep(:,merged_ROIs{i}(j)),d1,d2));
        title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
    end
    subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
    title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight;
    subplot(1,ln+2,ln+2);
    plot(1:T,(diag(max(C_keep(merged_ROIs{i},:),[],2))\C_keep(merged_ROIs{i},:))');
    hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
    title('Temporal Components','fontsize',16,'fontweight','bold')
    drawnow;
end

%% refine estimates excluding rejected components

Pm.p = p;    % restore AR value
[A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);


%% do some plotting

[A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
K_m = size(C_or,1);
Yr_ = single(reshape(im_stck,d,T)); % modified!!! 
[C_df,~] = extract_DF_F(Yr_,A_or,C_or,P_or,options); % extract DF/F values (optional) % modified!!! 

figure;
[Coor,json_file] = plot_contours(A_or,Cn,options,1); % contour plot of spatial footprints
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

%% display components

figure;
plot_components_GUI(Yr_,A_or,C_or,b2,f2,Cn,options);

%%
figure; 
imagesc(C_df)

%% WS Cleanup

if (clr_ws)
   clearvars A A2 A_keep Ain Am ans  b  bg_cmpt bidir_us bin C C2 ca_prfl C_keep ...
       center Cin Cm col_shift Coor clr_ws crop_image crop_coordinate d d1 d2 display_merging ...
       do_mtn_crcn do_scn_phs_crcn f fin fitness frm_rt gssn_hstm_stch_tshd ...
       gssn_xyz h hstm_stch_tshd i I in ind_cnn ind_corr ind_exc j K ...
       K_m keep M M1 manual_refine mk_mv merged_ROIs mn mx nb nFrames ...
       n_avrg n_frm n_nrn n_sgmt nrn_half options_rigid ...
       p P P2 Pm P_or phs_crcn_frm refine_components row1 row2 rval_space S S2...
       S_or shifts1 Sm T trgt_sgmt tau template1 th Y Yr YrA YrA2 YY im_stck
end 
%% Save workspace

if save_mat
    if brk_stck
        [path, filename, ~] = fileparts(tp_filename);
        fullfilename = path + "\" + filename + "-" + num2str(trgt_sgmt) + ".mat";
        save(fullfilename, '*', '-v7.3')
    else
        [path, filename, ~] = fileparts(tp_filename);
        fullfilename = path + "\" + filename + ".mat";
        save(fullfilename, '*', '-v7.3')
    end
end 
clear save_mat brk_stck
