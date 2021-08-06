function landwBud = readC2Vsim_LandwBud(filename, NsubRegions, Ntimes)
if isempty(NsubRegions)
    NsubRegions = 21;
end

if isempty(Ntimes)
    Ntimes = 1056;
end

Ag_Header = {'Area (AC)', 'Potential CUAW', 'Agricultural Supply Requirement',...
    'Pumping (-)', 'Diversion (-)', 'Shortage (=)', 'Re-use'}';
Urb_Header = Ag_Header;
Urb_Header{2} = 'Urban Supply Requirement';
Urb_Header(3,:) = [];
RIE_Header = {'Import', 'Export'};

frmt = '%s';
for i = 1:(length(Ag_Header) + length(Urb_Header) + length(RIE_Header))
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
		landwBud(isub,1).Time{i,1} = C{1,1}{1}(1:end-6);
        landwBud(isub,1).Ag.Header = Ag_Header;
        for j = 2:(length(Ag_Header)+1)
            landwBud(isub,1).Ag.Data(i,j-1) = C{1,j};
        end
        offset = (length(Ag_Header)+1);
        landwBud(isub,1).Urb.Header = Urb_Header;
        for j = 1:length(Urb_Header)
            landwBud(isub,1).Urb.Data(i,j) = C{1,j + offset};
        end
        offset = (length(Ag_Header) + length(Urb_Header) +1);
        landwBud(isub,1).RIE.Header = RIE_Header;
        for j = 1:length(RIE_Header)
            landwBud(isub,1).RIE.Data(i,j) = C{1,j + offset};
        end
    end
	temp = fgetl(fid);
    if exist ('OCTAVE_VERSION', 'builtin') == 0
        temp = fgetl(fid);
    end
end
fclose(fid);