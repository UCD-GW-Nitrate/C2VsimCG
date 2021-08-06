function streamBud = readC2Vsim_StreamBud(filename, NsubRegions, Ntimes)

if isempty(NsubRegions)
    NsubRegions = 21;
end

if isempty(Ntimes)
    Ntimes = 1056;
end

DataHeader = {'Upstream Inflow (+)', 'Downstream Outflow (-)', ...
    'Tributary Inflow (+)', 'Tile Drain (+)', 'Runoff (+)', ...
    'Return Flow (+)', 'Gain from Groundwater', 'Gain from Lake (+)', ...
    'Diversion (-)', 'By-pass Flow (-)','Discrepancy', 'Diversion Shortage'}';
frmt = '%s';
for i = 1:length(DataHeader) %This is the number of columns after the dates
	frmt = [frmt ' %f'];
end

fid = fopen(filename,'r');
for isub = 1:NsubRegions
	% skip 8 lines
    for i = 1:8
		temp = fgetl(fid);
    end
	
    %display(['SUBREGION ' num2str(isub)]);
    %if exist ('OCTAVE_VERSION', 'builtin') > 0
    %    fflush(stdout);
    %end
    
	for i = 1:Ntimes
		C = textscan(fid, frmt,1);
		streamBud(isub,1).Time{i,1} = C{1,1}{1}(1:end-6);
		for j = 2:(length(DataHeader)+1)
			streamBud(isub,1).Data(i,j-1) = C{1,j};
		end
	end
	streamBud(isub,1).Header = DataHeader;
	temp = fgetl(fid);
    if exist ('OCTAVE_VERSION', 'builtin') == 0
        temp = fgetl(fid);
    end
end
fclose(fid);