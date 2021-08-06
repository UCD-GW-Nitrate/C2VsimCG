function streamRchBud = readC2Vsim_StreamRchBud(filename, NsubRegions, Ntimes)

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
    for i = 1:7
        temp = fgetl(fid);
    end
	
    %display(['REACH ' num2str(isub)]);
    %if exist ('OCTAVE_VERSION', 'builtin') > 0
    %    fflush(stdout);
    %end
    
	for i = 1:Ntimes
		C = textscan(fid, frmt,1);
		streamRchBud(isub,1).Time{i,1} = C{1,1}{1}(1:end-6);
		for j = 2:(length(DataHeader)+1)
			streamRchBud(isub,1).Data(i,j-1) = C{1,j};
		end
	end
	streamRchBud(isub,1).Header = DataHeader;
	temp = fgetl(fid);
    if exist ('OCTAVE_VERSION', 'builtin') == 0
        temp = fgetl(fid);
    end
end
fclose(fid);