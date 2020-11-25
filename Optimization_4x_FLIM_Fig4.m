%% this is a code to super resolve experimental SPAD data 4x using the CMOS intensity and then find out the fluorescence lifetime imaging by tail fitting. 

%%%  The code has following parts:

%%%  Stage 1: This is the preprocessing of the data, which involves the removal of hot pixels,slided window binning and selection of the ROI from the spad data and CMOS data 

%%%  Stage 2: reconstruction for high resolution datacube

%%%  Stage 3: least square fitting to evaluate the lifetimes at each pixels
%%% of the low res raw data and the high res reconstructed data

% I_resize is the CMOS intensity loaded from the matfile CMOS_Raw.mat
% Spad is the low res SPAD Data loaded from the matfile Spad_Raw.mat

% The code uses functions SlideWinSum (needed prior to fitting) and
% decaymodelSingle (exponential decay model for tail fitting)


close all
clear
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PRE-PROCESSING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% read the intensity image

load ('C:\Users\areeba fatima\Dropbox\PoL in progress Shared\Presentations\work in progress\Areeba\Clara_data\ClaraPaper_revision\Final\Fig.4\CMOS_Raw.mat')        % this loads the CMOS intensity data as a variable named I_resize

% create a mask for CMOS intensity to remove background

I_medfilt = medfilt2(I_resize);
mask0 = I_medfilt;
mask0(mask0<175) = 0;
mask0(mask0>=175) = 1;
I_resize_filt = I_resize.*mask0;

% select the data segment to be used for optimization
I_final = I_resize_filt(153:456,77:228);

figure(1)
imagesc(I_final);
colormap jet
colorbar
axis equal tight
title('CMOS intensity')

I_final = double(I_final);     % this is the final CMOS intensity to be used in the optimization algorithm


%% Read the SPAD data histograms

load('C:\Users\areeba fatima\Dropbox\PoL in progress Shared\Presentations\work in progress\Areeba\Clara_data\ClaraPaper_revision\Final\Fig.4\Spad_Raw') % this loads the SPAD data measurement. The data is loaded as a variable named Spad
Spad_Intensity = (sum(Spad,3));

%% construct a mask,taking the filtered intensity as the reference.Use this mask to get the data from SPAD data

Int_medfilt = medfilt2(Spad_Intensity);
mask = Int_medfilt;
mask(mask<=100) = 0;
mask(mask>=100) = 1;

Spad = double(Spad);
for i=1:size(Spad,3)
    Spad_filt(:,:,i) = Spad(:,:,i).*mask;
end
Spad_Int_filt = sum(Spad_filt,3);

%% temporal binning. This is to reduce the timeframes from 200 to 50

k = 1;
for t = 1:4:200
    temp = Spad_filt(:,:,t:t+3);
    Spad_bin(:,:,k) = sum(temp,3);
    k = k+1;
end

%% Apply the slided window binning to each frame of raw data

for t = 1:50
    Spad_bin2(:,:,t)=SlideWinSum(Spad_bin(:,:,t),3);
end

%% select the Spad segment
Spad_bin3 = Spad_bin2(39:114,20:57,:);

%% resize the cmos intensity to match that of spad.This is to construct a reference for mask. Use the mask on SPAD data

I_resize_lr = imresize(I_final,[76,38],'nearest');

mask2 = I_resize_lr;
mask2(mask2<=0) = 0;
mask2(mask2>0) = 1;
mask2 = double(mask2);

for i=1:size(Spad_bin3,3)
    Spad_final(:,:,i) = Spad_bin3(:,:,i).*mask2;
end

Spad_final_Int = sum(Spad_final,3);

figure(2)
imagesc(Spad_final_Int)
colormap jet
axis equal tight
title('SPAD Intensity')


%%  %%%%%%%%%%%%%%%%%%%%%%%%%%% HIGH RES  RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct the blurring matrix
% The code process the data in four patches.Hence a blurring matrix of
% size 11552 x 11552 is constructed, corresponding to a single high res patch of size 76 x 152
 
band = 11;
sigma = 0.7;
Nr = 76;  %size of the rows
Nc = 152; %size of the column
zr = [exp(-((0:band-1).^2)/(2*sigma^2)),zeros(1,Nr-band)];
Ar = toeplitz(zr);
Ar = sparse(Ar);
zc = [exp(-((0:band-1).^2)/(2*sigma^2)),zeros(1,Nc-band)];
Ac = toeplitz(zc);
Ac = sparse(Ac);
A_b = (1/(2*pi*sigma^2))*kron(Ar,Ac);

%% optimization for high res dcube recovery
% the reconstruction is done in loop wherein in each iteration one of the
% patches is processed.

tic
p = 0;
for h = 0:19:57
   
    h1 = h+1;
    LR_Spad = Spad_final(h1:h1+18,1:38,:);
    h2 = 4*h+1;
    HR_Int = I_final(h2:h2+75,1:152);
    
    % increase the contrast of the cmos and scale it
    
    Spad_seg = sum(LR_Spad,3);

    a1 = sum(sum(HR_Int));
    a2 = 133*266;
    b1 = sum(sum(Spad_seg));
    b2 = 19*38;
    a = a1/a2;
    b = b1/b2;
    s = b/a;
    HR_Int_sc = HR_Int.*s;


   % start the optimization
   
   cvx_begin
   cvx_solver mosek
   variable x(76,152,50)
   Ux = x(1:end,2:end,:)-x(1:end,1:end-1,:);
   Uy = x(2:end,1:end,:)-x(1:end-1,1:end,:);
   Objective = norm((vec(DownsampleBlur_2(x,A_b)-LR_Spad)))+norm(vec(TempoIntegration(x)-HR_Int_sc))+ 1e-7*norm(vec(vec(x)),1)+1e-5*norm([Ux(:);Uy(:)],1)%+0.001*norm(SpatialInt(x)-16*K2)
   minimize(Objective)
   subject to 
    x>=0;
   cvx_end
   
   p = p+1;
   check{p} = x;
end
toc

%% combine the segments to get the full data

FullData = [check{1};check{2};check{3};check{4}];

% plot the intensity of the recovered high res data
INT = (sum(FullData,3));
figure(3)
imagesc(INT)
colormap jet
axis equal tight
colorbar
title('Reconstructed Data Intensity')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LIFETIME FITTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% least squares fitting for  the reconstructed data. The histograms are normalized to 1.Then the decay part of the histogram is fitted to a 
%  mono-exponential decay function (for which the function decaymodelSingle.m is required)



% define the time scale for the data

TBS_bin = (47*4)/1000;                % this is the timebin size in ns

tdata_full_t = TBS_bin.*(1:50);       % this is the time scale for the histogram (since there are 50 time bins)
fitRange2 = [15,50];                  % this is the range selected for tail fitting

tdata_bin = fitRange2(1):fitRange2(2);
tdata_bin = (TBS_bin.*tdata_bin);      %this is the time scale for the histogram (since there are 50 time bins)

initialcondition = [0,0,1,1];         % these are the intial parameters
lb = [0,0,0,0];                       % this is the lower bound
ub = [1,1,7,6];                       % this is the upper bound


% Least square fitting for all the pixels:

T_reco = zeros(304,152);           % this is to intialize a matrix that would store the lifetimes at all pixels of the reconstructed data

opts = optimset('Display','off');  % this prevents from the details to be printed for lsqcurvefit

for i = 1:304
    for j = 1:152
        if INT(i,j) > 4000        % this is the threshold with the intensity as reference
            ydata = squeeze(FullData(i,j,:));
            ydata_norm = ydata/max(ydata);
            ydata_n = (ydata_norm(15:50))';
            [xxx5,resnorm,residual,exitflag,output]=lsqcurvefit(@decaymodelSingle,initialcondition,tdata_bin,ydata_n,lb,ub,opts);
            T_reco(i,j) = xxx5(3);
        end
        
    end
    i
end

figure(4)
imagesc(T_reco)
colormap jet
axis equal tight
set(gca,'FontSize',14)
set(gca,'FontWeight','bold')
colorbar
caxis([0 3.0])
title('Lifetime recovered data')

mean_lifetime_reco= mean(nonzeros(T_reco)); 
std_lifetime_reco = std(nonzeros(T_reco)); 
fprintf('tau_reco = %.2f \x00B1 %.2f ns\n', mean_lifetime_reco,std_lifetime_reco)


%% least squares fitting for  the raw data


T_raw = zeros(76,38);                % this is to intialize a matrix that would store the lifetimes at all pixels of the reconstructed data


for i = 1:76
    for j = 1:38
        if Spad_final_Int(i,j) > 1300    % this is the threshold with the intensity as reference        
            ydata = squeeze(Spad_final(i,j,:));
            ydata_norm = ydata/max(ydata);
            ydata_n = (ydata_norm(15:50))';
            [xxx5,resnorm,residual,exitflag,output]=lsqcurvefit(@decaymodelSingle,initialcondition,tdata_bin,ydata_n,lb,ub,opts);
            T_raw(i,j) = xxx5(3);
        end
        
    end
    i
end

figure(5)
imagesc(T_raw)
colormap jet
axis equal tight
set(gca,'FontSize',14)
set(gca,'FontWeight','bold')
colorbar
caxis([0 3.0])
title('Lifetime raw data')



mean_lifetime_raw= mean(nonzeros(T_raw)); 
std_lifetime_raw = std(nonzeros(T_raw)); 
fprintf('tau_raw = %.2f \x00B1 %.2f ns\n', mean_lifetime_raw,std_lifetime_raw)



%% Lifetime Histograms

[sp,ed] = histcounts(nonzeros(T_raw),linspace(0,4,50));
for ii = 1:length(ed)-1
    ed2(ii) = (ed(ii)+ed(ii+1))/2;
end
[rc,ed] = histcounts(nonzeros(T_reco),ed);
figure(6)
bar(ed2,rc/max(rc))
hold on
bar(ed2,sp/max(sp))
legend('reconstruction','raw data')
xlim([1 3])
hold off


