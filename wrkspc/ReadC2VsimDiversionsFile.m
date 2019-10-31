function Diversions = ReadC2VsimDiversionsFile(filename)
fid = fopen(filename,'r');
section = 1;
cnt_steps = 1;
while ~feof(fid)
    tline = fgetl(fid);
    if strcmp(tline(1),'C')
        continue;
    end
    
    if section == 1
        c = textscan(tline,'%f');
        Ncoldiv = c{1};
        Diversions.header{1,1} = tline;
        for ii = 1:4
            tline = fgetl(fid);
            Diversions.header{1+ii,1} = tline;
        end
        section = section + 1;
        
    elseif section == 2
        c = textscan(tline,'%f/%f/%f');
        Diversions.Date(cnt_steps,1).Y = c{3};
        Diversions.Date(cnt_steps,1).M = c{2};
        Diversions.Date(cnt_steps,1).D = c{1};
        c = textscan(tline(18:end),'%f');
        Diversions.Data(cnt_steps,:) = c{1,1}';
        cnt_steps = cnt_steps + 1;
    end
    
end

fclose(fid);