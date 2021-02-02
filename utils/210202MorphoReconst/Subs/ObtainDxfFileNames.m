%%%
%%% Obtain 'dxf' filenames of a target directory.
%%% SortedFileNames are the output.
%%%

%%
%%
function SortedFileNames = ObtainDxfFileNames(FileDir)
%%
%%
	tmp         = dir(sprintf('%s\\*.dxf',FileDir));
	FileNames   = {tmp.name};
	FileNo      = [];
	for i = 1:numel(FileNames);
		DecompFile = textscan(FileNames{i},'%s %d %s', 'Delimiter', '.');
		FileNo = [FileNo, DecompFile{end-1}];
	end;
	[FileNo, I] = sort(FileNo);
	SortedFileNames = FileNames(I);

