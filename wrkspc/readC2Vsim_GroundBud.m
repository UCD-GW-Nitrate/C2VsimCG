function groundBud = readC2Vsim_GroundBud( filename, NsubRegions, Ntimes )
%readC2Vsim_GroundBud Read the CVground.BUD file
%   The script was written for a given version of IWFM and C2VSIM. It is
%   likely that other version would need modifications.
%   NsubRegions is the number of subregions. It sets to 21 if this is empty
%   Ntimes the number of timesteps to read. Default number is 1056

if isempty(NsubRegions)
    NsubRegions = 21;
end

if isempty(Ntimes)
    Ntimes = 1056;
end

fid = fopen(filename,'r');

frmt = '%s';
for i = 1:15 %This is the number of columns after the dates
	frmt = [frmt ' %f'];
end

DataHeader = {'Deep Percolation', 'Beginning Storage', 'Ending Storage', ...
	'Net Deep Percolation', 'Gain from Stream', 'Recharge', 'Gain from Lake', ...
	'Boundary Inflow', 'Subsidence', 'Subsurface Irrigation', 'Tile Drain Outflow', ...
	'Pumping', 'Net Subsurface Inflow', 'Discrepancy', 'Cumulative Subsidence'};

for isub = 1:NsubRegions
	% skip 8 lines
	for i = 1:8
		temp = fgetl(fid);
    end
	
     display(['Subregion ' num2str(isub)]);
    if exist ('OCTAVE_VERSION', 'builtin') > 0
        fflush(stdout);
    end
    
	for i = 1:Ntimes
		C = textscan(fid, frmt,1);
		groundBud(isub,1).Time{i,1} = C{1,1}{1}(1:end-6);
		for j = 2:16
			groundBud(isub,1).Data(i,j-1) = C{1,j};
		end
	end
	groundBud(isub,1).Header = DataHeader;
	temp = fgetl(fid);
    if exist ('OCTAVE_VERSION', 'builtin') == 0
        temp = fgetl(fid);
    end
end


fclose(fid);

end

