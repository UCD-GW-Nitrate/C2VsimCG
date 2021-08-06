function pS = runParetoSolutionsScen2(pS)

    for ii = 1:size(pS.x,1)
        display(['Solution ' num2str(ii)]);
        pS.Post(ii,1) = runSolution(pS.x(ii,:));
    end

end


function out = runSolution(x)
    % read existing diversion data 
    load('c2vsimDiversionData.mat','DivSpec','Diversions');
    
    % DATA Releated to the optimization input files
    % read optimization diversion time series
    DTS = readDTSfile(fullfile('..','OptimResults','inputFilesV2','DTS_APR21_Scen2_95_100.dat'));

    % read the economic function data
    elem_cost_data = readElemCostFile(fullfile('..','OptimResults','inputFilesV2','ElemCost_APR21_scen2.dat'));
    elem_cost_data.X_lift(elem_cost_data.X_lift > 0) = 0;
    
    % read discount rates
    disc_rates = readtable(fullfile('..','OptimResults','inputFilesV2','DiscountFactors3.dat'));
    
    %find the rivers that we divert from and append  the diversion data
    div_nodes = unique(elem_cost_data.DivNd(x == 1));
    nDivs = length(DivSpec);
    tot_land = 0;
    tot_capital = 0;
    tot_lift = 0;
    tot_conv = 0;
    for ii = 1:length(div_nodes)
        % Find the elements that receive diversion from this river node
        idts = find([DTS.Node]' == div_nodes(ii));
        Diversions.Data = [Diversions.Data DTS(idts,1).Ts];
        idx_el = find(elem_cost_data.DivNd == div_nodes(ii) & x' == 1);
        tmp_elem = elem_cost_data.ElemId(idx_el);
        DivSpec(nDivs + ii).ID = nDivs + ii;
        DivSpec(nDivs + ii).DESCR = ['Optim' num2str(ii)];
        DivSpec(nDivs + ii).IRDV = div_nodes(ii);
        DivSpec(nDivs + ii).ICDVMAX = 265;
        DivSpec(nDivs + ii).FDVMAX = 1;
        DivSpec(nDivs + ii).ICOLRL = size(Diversions.Data,2);
        DivSpec(nDivs + ii).FRACRL = 0.95;
        DivSpec(nDivs + ii).ICOLNL = size(Diversions.Data,2);
        DivSpec(nDivs + ii).FRACNL = 0.05;
        DivSpec(nDivs + ii).NDLDV = 1;
        DivSpec(nDivs + ii).IRGDL = 1;
        DivSpec(nDivs + ii).ICOLDL = size(Diversions.Data,2);
        DivSpec(nDivs + ii).FRACDL = 0;
        DivSpec(nDivs + ii).ICFSIRIG = 1;
        DivSpec(nDivs + ii).ICADJ = 1;
        DivSpec(nDivs + ii).NERELS = length(tmp_elem);
        DivSpec(nDivs + ii).IERELS = tmp_elem;
        DivSpec(nDivs + ii).FERELS = ones(length(tmp_elem),1);
        
        % Calculate the cost function
        tmpTS = DTS(idts,1).Ts(529:end);
        Qmax = max(tmpTS);
        totQ = sum(tmpTS);
        QmaxMAT(ii,1) = Qmax;
        QtotMAT(ii,1) = totQ;
        NreceivElem = length(idx_el);
        for j = 1:NreceivElem
            land = elem_cost_data.P_land(idx_el(j)) * ((Qmax/NreceivElem)*1000/75);
            tot_land = tot_land + land;
            capital = 5000*((Qmax/NreceivElem)*1000/75);
            tot_capital = tot_capital + capital;
            count_months = 1;
            cnt = 1;
            lift = 0;
            conveyance = 0;
            for k = 1:length(tmpTS)
                lift = lift + disc_rates.Var1(cnt) * ...
                    (0.17 * 1.45 * abs(elem_cost_data.X_lift(idx_el(j))) * (tmpTS(k)/NreceivElem) * 1000);
                conveyance = conveyance + disc_rates.Var1(cnt) * ...
                    (0.02 * elem_cost_data.X_dst(idx_el(j)) * (tmpTS(k)/NreceivElem) * 1000);
                count_months = count_months + 1;
                if count_months >= 13
                    count_months = 1;
                    cnt = cnt + 1;
                end
            end
            tot_lift = tot_lift + lift;
            tot_conv = tot_conv + conveyance;
        end
        
    end
    scenario_cost = tot_land + tot_capital + tot_lift + tot_conv;
    
    writeDiversionFiles(DivSpec, Diversions);
    currentdir = pwd;
    cd (fullfile('..','RunC2Vsim'));
    system([fullfile('..','c2vsim_cg_1921ic_r374_rev','C2VSim_CG_1921IC_R374_rev','bin','Simulation3.02.exe') ...
        ' CVsimPrintAll.in']);
    
    system([fullfile('..','c2vsim_cg_1921ic_r374_rev','C2VSim_CG_1921IC_R374_rev','bin','Budget3.02.exe') ...
        ' CVBudget.in']);
    cd(currentdir);
    Nsubregions = 21;
    Ntimes = 528;
    out.GWbud = readC2Vsim_GroundBud(fullfile('..','RunC2Vsim','Results','CVground.BUD'), Nsubregions+1, Ntimes);
    out.divbud = readC2Vsim_DivertBud(fullfile('..','RunC2Vsim','Results','CVdiverdtl.BUD'), Nsubregions, Ntimes);
    out.streamRchBud = readC2Vsim_StreamRchBud(fullfile('..','RunC2Vsim','Results','CVstreamrch.BUD'), 75, Ntimes);
    out.streamBud = readC2Vsim_StreamBud(fullfile('..','RunC2Vsim','Results','CVstream.BUD'), Nsubregions+1, Ntimes);
    out.rootBud = readC2Vsim_RootBud(fullfile('..','RunC2Vsim','Results','CVrootzn.BUD'), Nsubregions+1, Ntimes);
    out.landwBud = readC2Vsim_LandwBud(fullfile('..','RunC2Vsim','Results','CVlandwater.BUD'), Nsubregions+1, Ntimes);
    load('BasicScenario_DEC20.mat');
    out.EnvOF = out.GWbud(22).Data(end,3) - BasicBud.GWbud(22).Data(end,3);
    out.CostOF = scenario_cost;
    out.Qmax = QmaxMAT;
    out.Qtot = QtotMAT;
end

function TT = readElemCostFile(filename)
    T = readtable(filename);
    TT = renamevars(T,{'Var1','Var2','Var3','Var4','Var5','Var6','Var7'}, ...
        {'Id','DivNd','ElemId','Area','P_land','X_lift','X_dst'});
end

function DTS = readDTSfile(filename)
    fid = fopen(filename);
    CC = textscan(fid,'%f',2);
    Ndts = CC{1}(1);
    Ntsteps = CC{1}(2);
    for ii = 1:Ndts
        CC = textscan(fid,'%f',Ntsteps + 1);
        DTS(ii,1).Node = CC{1}(1);
        DTS(ii,1).Ts = CC{1}(2:end);
    end
    fclose(fid);
end

function writeDiversionFiles(DivSpec, Diversions)
    % Write diversion specification file
    fid = fopen(fullfile('..','RunC2Vsim','tempSpec.dat'),'w');
    fprintf(fid,'%d\n', length(DivSpec));
    for ii = 1:length(DivSpec)
        fprintf(fid, '%d %d %d %d %d %g %d %g %d %d %d %g %d %d\n', ...
            [DivSpec(ii,1).ID DivSpec(ii,1).IRDV DivSpec(ii,1).ICDVMAX DivSpec(ii,1).FDVMAX ...
            DivSpec(ii,1).ICOLRL DivSpec(ii,1).FRACRL DivSpec(ii,1).ICOLNL DivSpec(ii,1).FRACNL ...
            DivSpec(ii,1).NDLDV DivSpec(ii,1).IRGDL DivSpec(ii,1).ICOLDL DivSpec(ii,1).FRACDL ...
            DivSpec(ii,1).ICFSIRIG DivSpec(ii,1).ICADJ]);
    end
    
    for ii = 1:length(DivSpec)
        fprintf(fid, '%d %d %d %g\n', [DivSpec(ii,1).ID DivSpec(ii,1).NERELS DivSpec(ii,1).IERELS(1) DivSpec(ii,1).FERELS(1)]);
        if DivSpec(ii,1).NERELS > 1
            for j = 2:DivSpec(ii,1).NERELS
                fprintf(fid, '%d %g\n', [DivSpec(ii,1).IERELS(j) DivSpec(ii,1).FERELS(j)]);
            end
        end
    end
    fprintf(fid, '12\n60\n1min\n60\n1min\nC\n');
    load('c2vsimDiversionData.mat','BypassSpec');
    for ii = 1:length(BypassSpec.Data)
        fprintf(fid, '%d %d %d %d %g %g\n', ...
            [BypassSpec.Data(ii,1).ID, BypassSpec.Data(ii,1).IA BypassSpec.Data(ii,1).IDIVT ...
            BypassSpec.Data(ii,1).IDIVC BypassSpec.Data(ii,1).DIVRL BypassSpec.Data(ii,1).DIVNL]);
        if BypassSpec.Data(ii,1).IDIVC < 0
            fprintf(fid, '%f %f\n', [BypassSpec.Data(ii,1).DIVX BypassSpec.Data(ii,1).DIVY]');
        end
    end
    for ii = 1:length(BypassSpec.Data)
        fprintf(fid,'%d %d %d %g\n',[BypassSpec.Data(ii,1).ID length(BypassSpec.Data(ii,1).IERELS), ...
            BypassSpec.Data(ii,1).IERELS(1) BypassSpec.Data(ii,1).FERELS(1)]);
        if length(BypassSpec.Data(ii,1).IERELS) > 1
            fprintf(fid, '%d %g\n', [BypassSpec.Data(ii,1).IERELS(2:end) BypassSpec.Data(ii,1).FERELS(2:end)]');
        end
    end
    fclose(fid);
    
    
    % Write diversion data file
    fid = fopen(fullfile('..','RunC2Vsim','tempData.dat'), 'w');
    fprintf(fid,'%d\n', size(Diversions.Data,2));
    fprintf(fid,'%g\n', 43560000);
    fprintf(fid,'%d\n', 1);
    fprintf(fid,'%d\n', 0);
    fprintf(fid,'\n');
    for ii = 1:size(Diversions.Data,1)
        fprintf(fid, '%02g/%02g/%g_24:00', [Diversions.Date(ii,1).D Diversions.Date(ii,1).M Diversions.Date(ii,1).Y]);
        fprintf(fid, ' %g', Diversions.Data(ii,:));
        fprintf(fid, '\n');
    end
    fclose(fid);
end