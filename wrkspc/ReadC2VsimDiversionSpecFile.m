function [DivSpec, BypassSpec] = ReadC2VsimDiversionSpecFile(filename)
fid = fopen(filename,'r');


section = 1;
cnt_divs = 1;
tline_1 = [];
tline_2 = [];
while ~feof(fid)
    tline = fgetl(fid);
    if strcmp(tline(1),'C')
        tline_2 = tline_1;
        tline_1 = tline;
        continue;
    end
    
    if section == 1
        % Get the number of diversions
        c = textscan(tline,'%f');
        Ndivs = c{1};
        section = section + 1;
    elseif section == 2
        % get the list of diversions
        c = textscan(tline,'%f');
        DivSpec(cnt_divs, 1).ID = c{1}(1);
        DivSpec(cnt_divs, 1).DESCR = strtrim(tline_2(4:end));
        DivSpec(cnt_divs, 1).IRDV = c{1}(2);
        DivSpec(cnt_divs, 1).ICDVMAX = c{1}(3);
        DivSpec(cnt_divs, 1).FDVMAX = c{1}(4);
        DivSpec(cnt_divs, 1).ICOLRL = c{1}(5);
        DivSpec(cnt_divs, 1).FRACRL = c{1}(6);
        DivSpec(cnt_divs, 1).ICOLNL = c{1}(7);
        DivSpec(cnt_divs, 1).FRACNL = c{1}(8);
        DivSpec(cnt_divs, 1).NDLDV = c{1}(9);
        DivSpec(cnt_divs, 1).IRGDL = c{1}(10);
        DivSpec(cnt_divs, 1).ICOLDL = c{1}(11);
        DivSpec(cnt_divs, 1).FRACDL = c{1}(12);
        DivSpec(cnt_divs, 1).ICFSIRIG = c{1}(13);
        DivSpec(cnt_divs, 1).ICADJ = c{1}(14);
        cnt_divs = cnt_divs + 1;
        if cnt_divs > Ndivs
            section = section + 1;
            cnt_divs = 1;
        end
    elseif section == 3
        % Recharge zone for each diversion point
        c = textscan(tline,'%f');
        id = c{1}(1);
        Nel = c{1}(2);
        DivSpec(id, 1).NERELS = Nel;
        DivSpec(id, 1).IERELS = c{1}(3);
        DivSpec(id, 1).FERELS = c{1}(4);
        for ii = 2:Nel
            tline = fgetl(fid);
            c = textscan(tline,'%f');
            DivSpec(id, 1).IERELS = [DivSpec(id, 1).IERELS; c{1}(1)];
            DivSpec(id, 1).FERELS = [DivSpec(id, 1).FERELS; c{1}(2)];
        end
        cnt_divs = cnt_divs + 1;
        if cnt_divs > Ndivs
            section = section + 1;
            cnt_divs = 1;
        end
    elseif section == 4
        % Bypass Configuration Specifications
        c = textscan(tline,'%f');
        Nbypass = c{1};
        BypassSpec.header{1,1} = tline;
        for ii = 1:4
            tline = fgetl(fid);
            BypassSpec.header{1+ii,1} = tline;
        end
        section = section + 1;
    elseif section == 5
        c = textscan(tline,'%f');
        BypassSpec.Data(cnt_divs,1).ID = c{1}(1);
        BypassSpec.Data(cnt_divs,1).IA = c{1}(2);
        BypassSpec.Data(cnt_divs,1).IDIVT = c{1}(3);
        BypassSpec.Data(cnt_divs,1).IDIVC = c{1}(4);
        BypassSpec.Data(cnt_divs,1).DIVRL = c{1}(5);
        BypassSpec.Data(cnt_divs,1).DIVNL = c{1}(6);
        BypassSpec.Data(cnt_divs,1).DIVX = [];
        BypassSpec.Data(cnt_divs,1).DIVY = [];
        if BypassSpec.Data(cnt_divs,1).IDIVC < 0
            for ii = 1:abs(BypassSpec.Data(cnt_divs,1).IDIVC)
                tline = fgetl(fid);
                c = textscan(tline,'%f');
                BypassSpec.Data(cnt_divs,1).DIVX = [BypassSpec.Data(cnt_divs,1).DIVX; c{1}(1)];
                BypassSpec.Data(cnt_divs,1).DIVY = [BypassSpec.Data(cnt_divs,1).DIVY; c{1}(2)];
            end
        end
        cnt_divs = cnt_divs + 1;
        if cnt_divs > Nbypass
            section = section + 1;
            cnt_divs = 1;
        end
    elseif section == 6
        % Seepage locations for bypass canals
        c = textscan(tline,'%f');
        id = c{1}(1);
        Nel = c{1}(2);
        BypassSpec.Data(id,1).IERELS = c{1}(3);
        BypassSpec.Data(id,1).FERELS = c{1}(4);
        for ii = 2:Nel
            tline = fgetl(fid);
            c = textscan(tline,'%f');
            BypassSpec.Data(id,1).IERELS = [BypassSpec.Data(id,1).IERELS; c{1}(1)];
            BypassSpec.Data(id,1).FERELS = [BypassSpec.Data(id,1).FERELS; c{1}(2)];
        end
    end
    tline_2 = tline_1;
    tline_1 = tline;
end

fclose(fid);