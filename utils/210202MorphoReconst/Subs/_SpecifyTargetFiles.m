%%%
%%%
function SortedTargetFileNames = ...
	SpecifyTargetFiles(FileDir, SortedFileNames,TargetDomain)
%%%
%%%
	TargetFileNameID = [];
	for i = 1:numel(SortedFileNames);
%%%
%%%
		[namepoly, polyl, polynum] = ...
			DxfDataLoader(FileDir, char(SortedFileNames(i)));
		t = 0;
		for j = 1:polynum;
			for k = 1:numel(TargetDomain);
				if strcmp(char(namepoly{j}), char(TargetDomain{k}))
					t = t + 1;
				end;
			end;
		end;
		if (t > 0)
			TargetFileNameID = [TargetFileNameID,i];
		end;
%%%
%%%
	end;
%%%
%%%
	mmax = max(TargetFileNameID);
	mmin = min(TargetFileNameID);
	SortedTargetFileNames = SortedFileNames(mmin:mmax);
%%%
%%%

