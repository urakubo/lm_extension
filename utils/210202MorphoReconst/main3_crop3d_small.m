%%
%%
%%
clear;
targ = 5

addpath('./Subs');
p = ParamClass;
p.Addpaths;
p.SetTargetBranch(targ);
FILENAME = sprintf('%s%svoxels_rotated.mat', p.OutputDir, p.F);
load(FILENAME); % 'bw_Dend','bw_Mito','bw_PSD','bw_ER'

znum =  numel( bw_Dend(1,1,:) ) ;

max_row = zeros(znum,1)/0;
max_col = zeros(znum,1)/0;
min_row = zeros(znum,1)/0;
min_col = zeros(znum,1)/0;

for i = 1 : numel( bw_Dend(1,1,:) ) ;
	[row,col] = find(squeeze(bw_Dend(:,:,i)));
	if	length(row) == 0
		disp(i)
		continue;
	end
	max_row(i,1) = max(row);
	max_col(i,1) = max(col);
	min_row(i,1) = min(row);
	min_col(i,1) = min(col);
end

max_row = max(max_row);
max_col = max(max_col)-34;
min_row = min(min_row);
min_col = min(min_col);

bw_Dend_crop = bw_Dend([min_row:max_row],[min_col:max_col],:);
bw_Mito_crop = bw_Mito([min_row:max_row],[min_col:max_col],:);
bw_PSD_crop  = bw_PSD([min_row:max_row],[min_col:max_col],:);
bw_ER_crop   = bw_ER([min_row:max_row],[min_col:max_col],:);


% bw_Dend_crop(:,:,[[1:8],[825:882]]) = [];
% bw_Mito_crop(:,:,[[1:8],[825:882]]) = [];
% bw_PSD_crop(:,:,[[1:8],[825:882]])  = [];

bw_Dend_crop(:,:,[[1:392],[621:1068]]) = [];
bw_Mito_crop(:,:,[[1:392],[621:1068]]) = [];
bw_PSD_crop(:,:,[[1:392],[621:1068]])  = [];
bw_ER_crop(:,:,[[1:392],[621:1068]])   = [];


bw_Dend_crop = bw_Dend_crop([257:512],:,:);
bw_Mito_crop = bw_Mito_crop([257:512],:,:);
bw_PSD_crop = bw_PSD_crop([257:512],:,:);
bw_ER_crop  = bw_ER_crop([257:512],:,:);


figure('Name','cropped');
subplot(2,2,1);
imshow(squeeze(sum(bw_Dend_crop,1)));
subplot(2,2,2);
imshow(squeeze(sum(bw_Dend_crop,2)));
subplot(2,2,3);
imshow(squeeze(sum(bw_Dend_crop,3)));


figure('Name','cropped');
subplot(2,2,1);
imshow(squeeze(sum(bw_ER_crop,1)));
subplot(2,2,2);
imshow(squeeze(sum(bw_ER_crop,2)));
subplot(2,2,3);
imshow(squeeze(sum(bw_ER_crop,3)));


save('voxels_rotated_cropped_small.mat','bw_Dend_crop','bw_Mito_crop','bw_PSD_crop','bw_ER_crop');

bw_Dend_and_PSD = (bw_Dend_crop & bw_PSD_crop);
bw_Dend_not_PSD = (bw_Dend_crop & not(bw_PSD_crop));

hdf5write('CA1_small2.h5','/dendrite',uint8(bw_Dend_crop),'/PSD',uint8(bw_PSD_crop),'/Mitochondrion',uint8(bw_Mito_crop),'/ER',uint8(bw_ER_crop), '/dendrite_and_PSD', uint8(bw_Dend_and_PSD),'/dendrite_not_PSD', uint8(bw_Dend_not_PSD));



