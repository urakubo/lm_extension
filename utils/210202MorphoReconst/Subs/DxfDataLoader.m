%%%
%%% DxfDataLoader
%%%

function [namepoly, polylines, polynum] = DxfDataLoader(FileDir, FileName)

	%% Data load
	fname = sprintf('%s\\%s',FileDir,FileName);
	% disp(fname)
	fid   = fopen(fname);
	C     = textscan(fid,'%s','Delimiter','\n');
	fclose(fid);
	C=C{1,1};

	%% ID & NUM of Polyline
	idpoly     = strcmp('POLYLINE',C);
	polynum    = sum(idpoly);
	idpoly     = find(idpoly == 1);
	namepoly   = C(idpoly+8);

	%% ID of end of Polyline
	idployend  = strcmp('SEQEND',C);
	idployend  = find(idployend == 1);

	%% ID of Vertices
	idvert     = strcmp('VERTEX',C);
	idvert     = find(idvert == 1);

	%% preallocate variable
	polylines  = cell(polynum,1);

	%%
	for i=1:polynum;
	%%
		clear id targid xpoly ypoly polyline;
    	id     = (idvert > idpoly(i)) .* (idvert < idployend(i));
		targid = idvert(find(id));
		xpoly  = str2double(C(targid+2));
		ypoly  = str2double(C(targid+4));
    	polyline     = [xpoly ypoly]; % create polylinesubset    
    	polylines{i} = polyline;
	%%
	end;
	%%

%%%
%%% DataLoader
%%%

