
%%%
%%% Obtain Max Min
%%%

function [minX, maxX, minY, maxY] = MinMax(FileDir, SortedFileNames)

	minX = [];
	maxX = [];
	minY = [];
	maxY = [];
	for i = 1:numel(SortedFileNames);
		% fprintf('%s\n',char(SortedFileNames(i)));
		[namepoly, polylines, polynum] = ...
			DxfDataLoader(FileDir, char(SortedFileNames(i)));
		for j = 1:polynum;
			tmpX = polylines{j}(:,1);
			tmpY = polylines{j}(:,2);
 			minX = min([minX; tmpX]);
 			maxX = max([maxX; tmpX]);
 			minY = min([minY; tmpY]);
 			maxY = max([maxY; tmpY]);
		end;
	end;

