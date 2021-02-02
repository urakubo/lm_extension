%
% AquireVolumeImageFromSortedDxfFiles
%
% Requirements:
% DxfDataLoader
%


function [TargetDomainI, ZZ] = ...
		AquireVolumeImageFromSortedDxfFiles(FileDir, SortedTargetFileNames,...
			Canps, minX, minY, xypitch, zmult,TargetDomain);
	%%%
	TargetDomainI = [];
	for j = 1:zmult*4;
		TargetDomainI = cat(3, TargetDomainI, Canps);
	end;
	%%%
	for i = 1:numel(SortedTargetFileNames);
	%%%
		%%% Load data
		[namepoly, polyl, polynum] = ...
			DxfDataLoader(FileDir, char(SortedTargetFileNames(i)));
		
		% fprintf('Layer %s\n', char(SortedTargetFileNames(i))    );
		% namepoly 
		
		%%% Fill holes and binarization
		Cyt   = Canps;
		for j = 1:polynum;
			for k = 1:numel(TargetDomain) ; %% PSD
				if strcmp(char(namepoly{j}), TargetDomain{k}) %% PSD
					tmp = roipoly(Cyt,(polyl{j}(:,1)-minX)/xypitch, (polyl{j}(:,2)-minY)/xypitch );
					Cyt = Cyt + tmp;
				end;
			end;
		end;
		%%% Resample
		Cyt = (Cyt > 0.5);
		for j = 1:zmult;
			TargetDomainI = cat(3, TargetDomainI, Cyt);
		end;
	%%%
	end;
	%%%
	for j = 1:zmult*4;
		TargetDomainI = cat(3, TargetDomainI, Canps);
	end;
	ZZ = (0: numel(TargetDomainI(1,1,:))-1)*xypitch;
	%%%

