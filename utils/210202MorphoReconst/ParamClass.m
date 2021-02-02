%%
%%
classdef ParamClass < handle
%%
    properties
%%
%% I/O directories
%%
		F = filesep;
		IDIR        = 'C:\Users\uraku\Desktop\LatticeMicrobes\Laxmi_Morpho';
		ITARG_DIRS  = py.dict(pyargs('5', 'CA1_Dendrite5_mod'));
		MorphDir
		OutputDir

		zpitch  = 0.04;
		xypitch = 0.02;
		zmult   = 2;  %% round(obj.zpitch/obj.xypitch);
		RadiusPSD = 1;
%%
%% Target spine properties
%%
		SpineIDs
		HeadNameOrg
		HeadNames
		NeckNameOrg
		DendName
		Others
		TargetDomainNames
		TargetDomainNames_Cytosol
		TargetDomainNames_Mito
		TargetDomainNames_PSD
		TargetDomainNames_SER

%%
%% Paths to external programs
%%
		PATH_SUBROUTINES = 'Subs';
		DIRS_SUBROUTINES = {'iso2mesh','STLRead','FastMarching_version3b', ...
			'FastMarching_version3b\functions','FastMarching_version3b\shortestpath',...
			'distancePointLine','affine3d'};
%%
	end % Properties
%%
%%
	methods
		function obj = SetTargetBranch(obj,i)

			%% Set I/O directory

			obj.MorphDir    = sprintf('%s%s%s', obj.IDIR, obj.F, obj.ITARG_DIRS{num2str(i)} );

			%% Doamin names
			switch i
				case 2
					obj.SpineIDs    = 1:24;
					obj.HeadNameOrg = 'd2sp%02ghead';
					obj.NeckNameOrg = 'd2sp%02gneck';
					obj.DendName    = 'Dendrite2';
					obj.Others = {};
				case 3
					obj.SpineIDs    = 1:18;
					obj.HeadNameOrg = 'd3sp%02ghead';
					obj.NeckNameOrg = 'd3sp%02gneck';
					obj.DendName    = 'Dendrite3';
					obj.Others = {};
				case 4
					obj.SpineIDs    = 1:55;
					obj.HeadNameOrg = 'd4sp%02ghead';
					obj.NeckNameOrg = 'd4sp%02gneck';
					obj.DendName    = 'Dendrite4-2';
					obj.Others = {};
				case 5
					obj.SpineIDs    = 1:89;
					obj.HeadNameOrg = 'd5sp%02ghead';
					obj.NeckNameOrg = 'd5sp%02gneck';
					obj.DendName    = 'Dendrite5';
					obj.Others = {};
					obj.TargetDomainNames_Mito = {'Mitochondria1','Mitochondria2','Mitochondria3','Mitochondria4'};
					obj.TargetDomainNames_PSD  = {'PSD'};
					obj.TargetDomainNames_SER  = {'ER'};
					obj.OutputDir	= 'DATA_CA5';
				case 6
					obj.SpineIDs    = 1:27;
					obj.HeadNameOrg = 'd6sp%02ghead';
					obj.NeckNameOrg = 'd6sp%02gneck';
					obj.DendName    = 'Dendrite6';
					obj.Others = {'d6sp09-sp11branched neck'};
				case 7
					obj.SpineIDs    = 1:34;
					obj.HeadNameOrg = 'd7sp%02ghead';
					obj.NeckNameOrg = 'd7sp%02gneck';
					obj.DendName    = 'Dendrite7';
					obj.Others = {'d7sp17-18branched neck','d7sp20-21neck branched'};
				case 8
					obj.SpineIDs    = 1:34;
					obj.HeadNameOrg = 'd8sp%02ghead';
					obj.NeckNameOrg = 'd8sp%02gneck';
					obj.DendName    = 'Dendrite8-1';
					obj.Others = {};
				case 9
					obj.SpineIDs    = 1:20;
					obj.HeadNameOrg = 'd9sp%02ghead';
					obj.NeckNameOrg = 'd9sp%02gneck';
					obj.DendName    = 'Dendrite9';
					obj.Others = {};
				case 10
					obj.SpineIDs    = 1:21;
					obj.HeadNameOrg = 'd10sp%02ghead';
					obj.NeckNameOrg = 'd10sp%02gneck';
					obj.DendName    = 'Dendrite10';
					obj.Others = {};
   				otherwise
        			disp('Error! Check programs!')
			end
		end
%%
%%
		function Addpaths(obj)
%			for i = 1:numel(obj.DIRS_SUBROUTINES);
%				addpath(strcat(obj.PATH_SUBROUTINES, obj.F, obj.DIRS_SUBROUTINES{i}));
%			end;
			fontName = 'Arial';
			set(groot,'defaultAxesXColor','k'); % factory is [0.15,0.15,0.15]
			set(groot,'defaultAxesYColor','k');
			set(groot,'defaultAxesZColor','k');
			set(groot,'defaultAxesTickDir','out');
			set(groot,'defaultAxesBox','off');
			set(groot,'defaultAxesFontName', fontName);
			set(groot,'defaultTextFontName', fontName);
			set(groot,'defaultLegendFontName', fontName);
			set(groot,'defaultAxesFontName',     fontName);
			set(groot,'defaultTextFontName',     fontName);
			set(groot,'defaultLegendFontName',   fontName);
			set(groot,'defaultColorbarFontName', fontName);
		end
%%
%%
		function obj = SetTargetDomainNames(obj)
			obj.TargetDomainNames  = {};
			for i = 1:numel(obj.SpineIDs);
				HeadName    = sprintf(obj.HeadNameOrg, obj.SpineIDs(i));
				NeckName    = sprintf(obj.NeckNameOrg, obj.SpineIDs(i));
				obj.TargetDomainNames{2*i-1} = HeadName;
				obj.TargetDomainNames{2*i}   = NeckName;
			end;
			obj.TargetDomainNames{end+1} = obj.DendName;
			obj.TargetDomainNames        = [obj.TargetDomainNames, obj.Others];
			% disp(obj.TargetDomainNames)
		end;
%%
%%
%%
		function obj = SetTargetDomainNames_Cytosol(obj)
			obj.TargetDomainNames_Cytosol  = {};
			obj.HeadNames                  = {};
			for i = 1:numel(obj.SpineIDs);
				HeadName    = sprintf(obj.HeadNameOrg, obj.SpineIDs(i));
				obj.HeadNames{i} = HeadName;
			end;
			for i = 1:numel(obj.SpineIDs);
				HeadName    = sprintf(obj.HeadNameOrg, obj.SpineIDs(i));
				NeckName    = sprintf(obj.NeckNameOrg, obj.SpineIDs(i));
				obj.TargetDomainNames_Cytosol{2*i-1} = HeadName;
				obj.TargetDomainNames_Cytosol{2*i}   = NeckName;
			end;
			obj.TargetDomainNames_Cytosol{end+1} = obj.DendName;
			obj.TargetDomainNames_Cytosol        = [obj.TargetDomainNames_Cytosol, obj.Others];
			% disp(obj.TargetDomainNames_Cytosol)
		end;
%%
%%
%%
	end % methods
end
%%
%%


