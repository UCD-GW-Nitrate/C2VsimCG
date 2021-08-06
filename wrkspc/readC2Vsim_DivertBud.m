function divbud = readC2Vsim_DivertBud(filename, NsubRegions, Ntimes)
    str = fileread(filename);
    lines = regexp(str, '\r\n|\r|\n', 'split')';

    isub_line = 2;
    idiv_line = 5;
    istr_line = 6;
    ipm_line = 7;
    idt_line = 8;
    for ii = 1:NsubRegions
        % Parse Subregion ID
        C = strsplit(lines{isub_line,1});
        iSub = str2double(C{11});

        % Parse diversion Nodes
        C = strsplit(lines{idiv_line,1});
        C(:,1:3) = [];
        Div_ids = cellfun(@str2double,C)';
        Div_ids(isnan(Div_ids),:) = [];

        % Parse stream Nodes
        C = strsplit(lines{istr_line,1});
        C(:,1:4) = [];
        Strm_ids = cellfun(@str2double,C)';
        Strm_ids(isnan(Strm_ids),:) = [];

        % Parse plus minus
        C = strsplit(lines{ipm_line,1});
        pm = cellfun(@isplus, C)';
        pm(isnan(pm),:) = [];
        data = nan(Ntimes, 2*length(Div_ids));
        YMD = nan(Ntimes,3);
        for j = 1:Ntimes
            C = strsplit(lines{idt_line+j,1})';
            cctime = textscan(C{1,1},'%f/%f/%f/_%s');
            C(1,:) = [];
            data(j,:) = cellfun(@removeBrakets, C)';
            YMD(j,:) = [cctime{3} cctime{1} cctime{2}];
        end
        divbud(iSub,1).ID = iSub;
        divbud(iSub,1).DivID = Div_ids;
        divbud(iSub,1).RivID = Strm_ids;
        divbud(iSub,1).Plus = pm;
        divbud(iSub,1).Data = data;
        divbud(iSub,1).YMD = YMD;
        isub_line = isub_line + 6 + Ntimes + 3;
        idiv_line = idiv_line + 3 +  Ntimes + 6;
        istr_line = istr_line + 2 +  Ntimes + 7;
        ipm_line = ipm_line + 1 + Ntimes + 8;
        idt_line = idt_line + Ntimes + 9;
    end

end

function p = isplus(c)
    if isempty(c)
        p = nan;
    else
        if strcmp(c(2),'+')
            p = 1;
        elseif strcmp(c(2),'-')
            p = -1;
        else
            p = nan;
        end
    end    
end
function c = removeBrakets(c)
    c = str2double(c(1:end-1));
end