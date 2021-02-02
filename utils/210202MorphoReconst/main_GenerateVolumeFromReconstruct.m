%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mod by H Urakubo                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%%
function main_GenerateVolumeFromReconstruct
%%%
%%%

	targ = 5
	savename = 'CA1_dend5'

	addpath('./Subs');
	p = ParamClass;
	p.Addpaths;
	p.SetTargetBranch(targ);
	p.SetTargetDomainNames_Cytosol;
	disp(p)
	mkdir(p.OutputDir);


%%%
%%% Obtain 'dxf' filenames in a target directory.
%%%

	SortedFileNames  = ObtainDxfFileNames(p.MorphDir);


%%%
%%% Create canpus
%%%

	fprintf('OK... (line 27) \n');
  	[minX, maxX, minY, maxY] = MinMax(p.MorphDir, SortedFileNames);
	fprintf('Min X: %f\n', minX);
	fprintf('Max X: %f\n', maxX);
	fprintf('Min Y: %f\n', minY);
	fprintf('Max Y: %f\n', maxY);

	XX = 0:p.xypitch:maxX;
	YY = 0:p.xypitch:maxY;

	tmp = zeros(numel(XX),numel(YY));
	size(tmp)
	Canps = mat2gray(tmp);

%%%
%%% Volume data aquition
%%%

%	p.TargetDomainNames_Cytosol

	[DendSpineI, ZZ] = ...
		AquireVolumeImageFromSortedDxfFiles(p.MorphDir, SortedFileNames, ...
			Canps, minX, minY, p.xypitch, p.zmult, p.TargetDomainNames_Cytosol);

	[mitochondriaI, ZZ] = ...
		AquireVolumeImageFromSortedDxfFiles(p.MorphDir, SortedFileNames, ...
			Canps, minX, minY, p.xypitch, p.zmult, p.TargetDomainNames_Mito);

	[ERI, ZZ] = ...
		AquireVolumeImageFromSortedDxfFiles(p.MorphDir, SortedFileNames, ...
			Canps, minX, minY, p.xypitch, p.zmult, p.TargetDomainNames_SER);

	[PSDI, ZZ] = ...
		AquireVolumeImageFromSortedDxfFiles(p.MorphDir, SortedFileNames, ...
			Canps, minX, minY, p.xypitch, p.zmult, p.TargetDomainNames_PSD);

	% DendSpineI = SimpleSmoothing(DendSpineI, p.Radius);

	StackedI = uint8((DendSpineI) > 0.5) + uint8(mitochondriaI + ERI > 0.5);
	bw_Dend  = logical(DendSpineI > 0.5);
	bw_Mito  = logical(mitochondriaI > 0.5);
	bw_ER    = logical(ERI > 0.5);
	bw_PSD   = logical(PSDI > 0.5);


%%%
%%% Smoothing & Dilution
%%%

	PSDI  = SimpleSmoothing(PSDI, p.RadiusPSD);


%%%
%%% Plot a target spine
%%%

%%{
	fv_Dend   = isosurface(YY,XX,ZZ,DendSpineI,0.5);
	fv_Mito   = isosurface(YY,XX,ZZ,bw_Mito,0.5);
	fv_PSD    = isosurface(YY,XX,ZZ,bw_PSD ,0.5);

	figure;
	set(gca,'DataAspectRatio',[1 1 1]);
	p1    = patch(fv_Mito,'FaceColor','g','EdgeColor','none','FaceAlpha',.5);
	hold on;
	p1    = patch(fv_Dend,'FaceColor','b','EdgeColor','none','FaceAlpha',.5);
	p1    = patch(fv_PSD ,'FaceColor','r','EdgeColor','none','FaceAlpha',.5);

%	FILENAME = sprintf('./%s/%s_1.fig',p.OutputDir,dd);
%	saveas(gcf,FILENAME);
%	return;
%%}



%%%
%%% Save Data
%%%


	FILENAME = sprintf('./%s/%s.mat',p.OutputDir,savename);
	save(FILENAME,'-v7.3');

	FILENAME = sprintf('%s%svoxels.mat', p.OutputDir, p.F)
	save(FILENAME,'bw_Dend','bw_Mito','bw_PSD','bw_ER');

%%%
%%% Smoothing
%%%

function StackedI = SimpleSmoothing(StackedI, Radius);

	Stacked1st = logical(StackedI > 0.5);
	Stacked1st = imdilate(Stacked1st, strel('sphere', Radius));
	Stacked1st = imerode(Stacked1st, strel('sphere', Radius));
	StackedI   = uint8(Stacked1st > 0.5);


%%%
%%% Save Elem
%%%

function SaveElem(FILENAME, elem_);
	fileID = fopen(FILENAME,'wt');
	fprintf(fileID,'%d %d %d %d\n', elem_');
	fclose(fileID);


%%%
%%% Save ElemID
%%%

function SaveElemID(FILENAME, elemID);
	fileID = fopen(FILENAME,'wt');
	fprintf(fileID,'%d\n',  elemID);
	fclose(fileID);


%%%
%%%
%%%


