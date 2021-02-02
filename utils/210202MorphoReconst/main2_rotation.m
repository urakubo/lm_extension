%%
%%
%%
%%
clear

targ = 5

addpath('./Subs');
p = ParamClass;
p.Addpaths;
p.SetTargetBranch(targ);
FILENAME = sprintf('%s%svoxels.mat', p.OutputDir, p.F)
load(FILENAME); % 'bw_Dend','bw_Mito','bw_PSD','bw_ER'

sizeIn = size(bw_Dend);
hFigOriginal = figure;
hAxOriginal  = axes;
slice(double(bw_Dend),sizeIn(2)/5,sizeIn(1)/5,sizeIn(3)/5);
grid on, shading interp, colormap gray


%{
Z axis rotation
theta = -pi/4;
t = [cos(theta) , -sin(theta),	0, 0;
	sin(theta)	, cos(theta) ,  0, 0;
	0			, 0		     ,	1, 0;
	0			, 0			 ,	0, 1];
%}

% Y axis rotation
theta = -pi/15;
tt = [cos(theta) , 0, -sin(theta),	0;
	0			, 1,	0		,	0;
	sin(theta)	, 0, cos(theta)	,	0;
	0			, 0,	0		,	1];

% X axis rotation
theta = -(pi-0.3)/3;
t = [	1		, 0,	0		,	0;
    	0, cos(theta) , -sin(theta),	0;
		0, sin(theta) , cos(theta)	,	0;
		0			, 0,	0		,	1];
t = tt*t;

tform = affine3d(t);
DendRotated = imwarp(double(bw_Dend),tform);
MitoRotated = imwarp(double(bw_Mito),tform);
PSDRotated  = imwarp(double(bw_PSD),tform);
ERRotated  = imwarp(double(bw_ER),tform);

%size(DendRotated)
%size(MitoRotated)
%size(PSDRotated)

% figure('Name','org');subplot(2,2,1);imshow(squeeze(sum(bw_Dend,1)));;subplot(2,2,2);imshow(squeeze(sum(bw_Dend,2)));subplot(2,2,3);imshow(squeeze(sum(bw_Dend,3)));

figure('Name','rotated');subplot(2,2,1);imshow(squeeze(sum(DendRotated,1)));;subplot(2,2,2);imshow(squeeze(sum(DendRotated,2)));subplot(2,2,3);imshow(squeeze(sum(DendRotated,3)))

figure('Name','ERRotated');subplot(2,2,1);imshow(squeeze(sum(ERRotated,1)));;subplot(2,2,2);imshow(squeeze(sum(ERRotated,2)));subplot(2,2,3);imshow(squeeze(sum(ERRotated,3)))


bw_Dend = logical(DendRotated);
bw_Mito = logical(MitoRotated);
bw_PSD  = logical(PSDRotated);
bw_ER   = logical(ERRotated);


FILENAME = sprintf('%s%svoxels_rotated.mat', p.OutputDir, p.F)
save(FILENAME,'bw_Dend','bw_Mito','bw_PSD','bw_ER');



