function rootBud = readC2Vsim_RootBud(filename, NsubRegions, Ntimes)

if isempty(NsubRegions)
    NsubRegions = 21;
end

if isempty(Ntimes)
    Ntimes = 1056;
end

Ag_UrbHeader = {'Area (AC)', 'Precipitation', 'Runoff',...
    'Prime Applied Water', 'Reused Water', 'Total Applied Water', ...
    'Return Flow', 'Beginning Storage', 'Net Gain from Land Expansion (+)', ...
    'Infiltration (+)', 'Actual ET (-)','Deep Percolation (-)', 'Ending Storage (=)'}';
NV_RP_Header = Ag_UrbHeader;
NV_RP_Header([4:7],:) = [];
frmt = '%s';
for i = 1:(2*length(Ag_UrbHeader) + length(NV_RP_Header)) %This is the number of columns after the dates
	frmt = [frmt ' %f'];
end

fid = fopen(filename,'r');
for isub = 1:NsubRegions
	% skip 8 lines
    for i = 1:10
        temp = fgetl(fid);
    end
    
    %display(['SUBREGION ' num2str(isub)]);
    %if exist ('OCTAVE_VERSION', 'builtin') > 0
    %    fflush(stdout);
    %end
    
    for i = 1:Ntimes
		C = textscan(fid, frmt,1);
		rootBud(isub,1).Time{i,1} = C{1,1}{1}(1:end-6);
        rootBud(isub,1).Ag.Header = Ag_UrbHeader;
        for j = 2:(length(Ag_UrbHeader)+1)
            rootBud(isub,1).Ag.Data(i,j-1) = C{1,j};
        end
        offset = (length(Ag_UrbHeader)+1);
        rootBud(isub,1).Urb.Header = Ag_UrbHeader;
        for j = 1:length(Ag_UrbHeader)
            rootBud(isub,1).Urb.Data(i,j) = C{1,j + offset};
        end
        offset = (2*length(Ag_UrbHeader)+1);
        rootBud(isub,1).NVRP.Header = Ag_UrbHeader;
        for j = 1:length(NV_RP_Header)
            rootBud(isub,1).NVRP.Data(i,j) = C{1,j + offset};
        end
    end
	temp = fgetl(fid);
    if exist ('OCTAVE_VERSION', 'builtin') == 0
        temp = fgetl(fid);
    end
end
fclose(fid);