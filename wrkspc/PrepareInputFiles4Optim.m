%% Model path
c2vsim_path = ['..' filesep 'c2vsim_cg_1921ic_r374_rev' filesep 'C2VSim_CG_1921IC_R374_rev' filesep];
load('diversionRiverNodes')
ScenID = 1;
scen = DiversionScenarios(['Scen' num2str(ScenID)]);
%% Remove Chowchilla and Fresno
% We remove these rivers because their stream hydrographs are negligible
% compared to the other selected rivers
scen([5 6],:) = [];
%% Find on how many elements each node can recharge
for ii = 1:length(scen)
   n = 0;
   scen(ii,1).IE = [];
   for j = 1:length(scen(ii,1).DivIds)
       scen(ii,1).IE = [scen(ii,1).IE;diversionRiverNodes(scen(ii,1).DivIds(j,1),1).IE];
   end
   scen(ii,1).IE = unique(scen(ii,1).IE);
   scen(ii,1).Ne = length(scen(ii,1).IE);
end
%% Make a unique list of elements
IE_unique = [];
for ii = 1:length(scen)
    IE_unique = [IE_unique; scen(ii,1).IE];
end
%% Find the elements that can receive water from multiple diversion nodes
sorted_IE = sort(IE_unique);
IE_multDiv = sorted_IE(diff(sorted_IE) == 0);
IE_multDiv = mat2cell(IE_multDiv, ones(1,length(IE_multDiv)),1);
for ii = 1:length(IE_multDiv)
    IE_multDiv{ii,2} = [];
    for j = 1:length(scen)
       id = find(scen(j,1).IE == IE_multDiv{ii,1});
       if ~isempty(id)
           IE_multDiv{ii,2} = [IE_multDiv{ii,2};j];
       end
    end
end
%% Manually delete the elements
% For Scen 1 we delete the nodes here and set manual_select_IE as empty
scen(6).IE([2 4 6 15 17],:) = [];
scen(5).IE([11 12 13],:) = [];
%% Manualy choose the diversion node for those elements
manual_select_IE = cell2mat(IE_multDiv(:,1));
manual_select = [8;8;8;8;4;5;7;7;7;7;7;7];
%% or set as empty
manual_select_IE = [];
manual_select = [];
%% Read node coordinates and elevations
c2vsim_nodes = shaperead(['..' filesep 'gis_data' filesep 'C2Vsim_nodes_3310']);
ND = [[c2vsim_nodes.ID]' [c2vsim_nodes.X]' [c2vsim_nodes.Y]' [c2vsim_nodes.ELV]'];
c2vsim_river_nodes = shaperead(['..' filesep 'gis_data' filesep 'C2Vsim_riverNodes.shp']);
c2vsim_rivers = shaperead(['..' filesep 'gis_data' filesep 'C2Vsim_rivers_3310.shp']);
IGW_IRV = [[c2vsim_river_nodes.IGW]' [c2vsim_river_nodes.IRV]'];
fid = fopen([c2vsim_path 'Preprocessor' filesep 'CVelement.dat'],'r');
temp = textscan(fid, '%f %f %f %f %f', 1392, 'HeaderLines',93);
fclose(fid);
MSH = [temp{1,2} temp{1,3} temp{1,4} temp{1,5}];
%%
% 1) RiverName : The name of the river where the water is diverted from
% 2)  DivNdId :  River node ID (IRV) where the water is extracted from the river in C2Vsim simulation
% 3)  DivNdElev : The elevation of the above node
% 4)  DivNdOpt : The river node id where the water is extracted from the river during optimization.
% 5)  ElemID: The element id that receives water from the diversion node
% 6)  ElemElevAv : The average element elevation (average of the element corners)
% 7)  ElemElevMax : The elevation of the element corner with maximum elevation
% 8)  ElemElevMin : The elevation of the element corner with minimum elevation
% 9)  DistAct: The actual distance between Elem barycenter and the DivNdId
% 10) DistNear: The distance between the ElemID barycenter and the closest stream node of the RiverName stream
% 11) NearNdElev : The elevation of the closest stream node of the river 
IE_unique = unique(IE_unique);
DATA = cell(0,11);
cnt = 1;
for ii = 1:length(IE_unique)
    tmp_data = cell(1,11);
    % check if it is one of the manually selected
    test = find(manual_select_IE == IE_unique(ii));
    if isempty(test)
        for j = 1:length(scen)
            test = find(scen(j,1).IE == IE_unique(ii));
            if ~isempty(test)
                break
            end
        end
    else
        j = manual_select(test);
    end
    tmp_data{1,1} = scen(j,1).RiverName;
    
    for k = 1:length(scen(j,1).DivIds)
        test = find(diversionRiverNodes(scen(j,1).DivIds(k),1).IE == IE_unique(ii));
        if ~isempty(test)
            if diversionRiverNodes(scen(j,1).DivIds(k),1).IRDV == 0
                tmp_data{1,2} = scen(j,1).ImportNode;
            else
                tmp_data{1,2} = diversionRiverNodes(scen(j,1).DivIds(k),1).IRDV;
            end
            break;
        end
    end
    if isempty(test)
        warning(['Cant find diversion node for element ' + num2str(IE_unique(ii))]);
    end
    
    
    igw = IGW_IRV(find(IGW_IRV(:,2) == tmp_data{1,2},1));
    tmp_data{1,3} = ND(igw,4);
    tmp_data{1,4} = scen(j,1).NodeInSim;
    tmp_data{1,5} = IE_unique(ii);
    el_nd = MSH(IE_unique(ii),:);
    el_nd(:,el_nd == 0) = [];
    tmp_data{1,6} = mean(ND(el_nd,4));
    tmp_data{1,7} = max(ND(el_nd,4));
    tmp_data{1,8} = min(ND(el_nd,4));
    el_barycent = mean(ND(el_nd,2:3));
    tmp_data{1,9} = sqrt(sum((el_barycent - ND(igw,2:3)).^2))/1000;
    
    
    min_dst = 99999999999999;
    for k = 1:length(scen(j,1).rivID)
        for kk = 1:length(c2vsim_rivers(scen(j,1).rivID(k),1).X) - 1
            tmp_dst = sqrt((el_barycent(1) - c2vsim_rivers(scen(j,1).rivID(k),1).X(kk))^2 + (el_barycent(2) - c2vsim_rivers(scen(j,1).rivID(k),1).Y(kk))^2);
            if tmp_dst < min_dst
                min_dst = tmp_dst;
                xmin = c2vsim_rivers(scen(j,1).rivID(k),1).X(kk);
                ymin = c2vsim_rivers(scen(j,1).rivID(k),1).Y(kk);
            end
        end
    end
    if ~isempty(scen(j,1).rivID)
        tmp_data{1,10} = min_dst/1000;
        [cc, dd] = min(sqrt((ND(:,2) - xmin).^2 + (ND(:,3) - ymin).^2));
        tmp_data{1,11} = ND(dd,4);
    end
    
    %elem_divND_dist = [];
    %for k = 1:length(el_nd)
    %    elem_divND_dist = [elem_divND_dist;
    %        sqrt(sum((ND(el_nd(k),2:3) - ND(igw,2:3)).^2)) ND(MSH(IE_unique(ii),k),4)];
    %end
    
    %[cc, dd] = min(elem_divND_dist(:,1));
    %tmp_data{1,10} = cc/1000;
    %tmp_data{1,11} = elem_divND_dist(dd,2);
    DATA = [DATA;tmp_data]; 
end
%% append Header and write
DATA_HEADER = {'RiverName','DivNdId', 'DivNdElev', 'DivNdOpt','ElemID','ElemElevAv', 'ElemElevMax', 'ElemElevMin', 'DistAct', 'DistNear', 'NearNdElev'};
writecell([DATA_HEADER;DATA], ['Scenario' num2str(ScenID) '_OCT20'],'FileType','spreadsheet');
%% Create shapefiles with the diversions data
% First a point shapefile with the diversion nodes
% and an polygon shapefile with attributes the node where the diversion
% takes place.
%% Node shapefile
clear S
S(1,1).Geometry = 'Point';
S(1,1).X = [];
S(1,1).Y = [];
S(1,1).IRV = [];
S(1,1).IGW = [];
S(1,1).RivName = [];
S(1,1).IRVOptm=[];
Unique_div_nodes = unique(cell2mat(DATA(:,2)));
for ii = 1:length(Unique_div_nodes)
    % find igw
    igw = IGW_IRV(find(IGW_IRV(:,2) == Unique_div_nodes(ii),1));
    
    
    S(ii,1).Geometry = 'Point';
    S(ii,1).X = ND(igw,2);
    S(ii,1).Y = ND(igw,3);
    S(ii,1).IRV = Unique_div_nodes(ii);
    S(ii,1).IGW = igw;
    test = find([scen.ImportNode]' == Unique_div_nodes(ii));
    if ~isempty(test)
        S(ii,1).RivName = scen(test,1).RiverName;
        S(ii,1).IRVOptm = scen(test,1).NodeInSim;
    else
        for j = 1:length(scen)
            for k = 1:length(scen(j,1).DivIds)
               if diversionRiverNodes(scen(j,1).DivIds(k),1).IRDV == Unique_div_nodes(ii)
                   S(ii,1).RivName = scen(j,1).RiverName;
                   S(ii,1).IRVOptm = scen(j,1).NodeInSim;
                   break;
               end
            end
        end
    end
end
shapewrite(S, ['..' filesep 'gis_data' filesep 'Scen' num2str(ScenID) 'DivNodes_OCT20_3310']);
%% Polygon shapefile
clear SS
SS(1,1).Geometry = 'Polygon';
SS(1,1).X = [];
SS(1,1).Y = [];
SS(1,1).BoundingBox = [];
SS(1,1).IE = [];
SS(1,1).DivND = [];
SS(1,1).DivNDOptim = [];
SS(1,1).RivName = [];

for ii = 1:size(DATA,1)
    SS(ii,1).Geometry = 'Polygon';
    el_nd = MSH(DATA{ii,5},:);
    el_nd(:,el_nd == 0) = [];
    el_nd = [el_nd el_nd(1)];
    SS(ii,1).X = [ND(el_nd,2)' nan];
    SS(ii,1).Y = [ND(el_nd,3)' nan];
    SS(ii,1).BoundingBox = [min(ND(el_nd,2)) min(ND(el_nd,3));max(ND(el_nd,2)) max(ND(el_nd,3))];
    SS(ii,1).IE = DATA{ii,5};
    SS(ii,1).DivND = DATA{ii,2};
    SS(ii,1).DivNDOptim = DATA{ii,4};
    SS(ii,1).RivName = DATA{ii,1};
end
shapewrite(SS, ['..' filesep 'gis_data' filesep 'Scen' num2str(ScenID) 'DivElem_OCT20_3310']);
%% Plot the stream hydrographs for the selected diversion nodes
load('c2vsimBasicScenario.mat', 'SWhyd')
%%
start_date = datetime(1965,10,1);
tm = start_date + calmonths(0:527);
colors = colororder('default');
fact = 1233.48*1e-9;
%%
figure(10);
clf
ii = 3;
this_stream_hyd = SWhyd.data(529:1056,scen(ii).NodeInSim);
prc_95 = prctile(this_stream_hyd,95);
prc_90 = prctile(this_stream_hyd,90);
plot(tm, fact * this_stream_hyd, 'DisplayName', 'Hydrograph','color',colors(1,:),'linewidth',2)
hold on
plot(tm([1,end]), fact*prc_95*ones(2,1),'linewidth',2,'DisplayName','95^{th} percentile','color',colors(2,:))
plot(tm([1,end]), fact*prc_90*ones(2,1),'linewidth',2,'DisplayName','90^{th} percentile','color',colors(3,:))
plot(tm([1,end]), fact*(prc_95 + 100000)*ones(2,1),'--','linewidth',2,'DisplayName','95^{th} 100 TAF CAP','color',colors(2,:))
plot(tm(1:40:length(tm)), fact*(prc_95 + 200000)*ones(length(1:40:length(tm)),1),':','linewidth',2,'DisplayName','95^{th} 200 TAF CAP','color',colors(2,:))
plot(tm([1,end]), fact*(prc_90 + 100000)*ones(2,1),'--','linewidth',2,'DisplayName','90^{th} 100 TAF CAP','color',colors(3,:))
plot(tm(1:40:length(tm)), fact*(prc_90 + 200000)*ones(length(1:40:length(tm)),1),':','linewidth',2,'DisplayName','90^{th} 200 TAF CAP','color',colors(3,:))

title(['Hydrograph for ' scen(ii).RiverName ' river'])
ylabel('Water volume [km^3]')
xlabel('Time')
legend('Location','northwest')
grid on
% print -dpng -r600 Tuolumne_hydrograph
%%
figure(20);
clf
hyd_95 = zeros(length(tm),1);
hyd_90 = zeros(length(tm),1);
hyd_95(this_stream_hyd > prc_95) = this_stream_hyd(this_stream_hyd > prc_95) - prc_95;
hyd_90(this_stream_hyd > prc_90) = this_stream_hyd(this_stream_hyd > prc_90) - prc_90;
plot(tm, fact*cumsum(min(hyd_95,100000)), 'DisplayName', '95^{th} 100 CAP','linewidth',2,'color',colors(2,:)) 
hold on
plot(tm, fact*cumsum(min(hyd_95,200000)),'--', 'DisplayName', '95^{th} 200 CAP','linewidth',2,'color',colors(2,:)) 
plot(tm, fact*cumsum(min(hyd_90,100000)), 'DisplayName', '90^{th} 100 CAP','linewidth',2,'color',colors(3,:)) 
plot(tm, fact*cumsum(min(hyd_90,200000)),'--', 'DisplayName', '90^{th} 200 CAP','linewidth',2,'color',colors(3,:))
title(['Cumulative diversion for ' scen(3).RiverName ' river'])
ylabel('Water volume [km^3]')
xlabel('Time')
legend('Location','northwest')
grid on
% print -dpng -r600 Tuolumne_cumulDiversion
%% plot Cumulative diversions for all scenarions in panels
cumDiv_95_100 = zeros(length(scen),length(tm));
cumDiv_95_200 = zeros(length(scen),length(tm));
cumDiv_90_100 = zeros(length(scen),length(tm));
cumDiv_90_200 = zeros(length(scen),length(tm));
for ii = 1:length(scen)
    str_hyd = SWhyd.data(529:1056,scen(ii).NodeInSim);
    prc_95 = prctile(str_hyd,95);
    prc_90 = prctile(str_hyd,90);
    
    tmp = zeros(length(tm),1);
    tmp(str_hyd > prc_95) = str_hyd(str_hyd > prc_95) - prc_95;
    cumDiv_95_100(ii,:) = min(tmp, 100000);
    cumDiv_95_200(ii,:) = min(tmp, 200000);
    
    tmp = zeros(length(tm),1);
    tmp(str_hyd > prc_90) = str_hyd(str_hyd > prc_90) - prc_90;
    cumDiv_90_100(ii,:) = min(tmp, 100000);
    cumDiv_90_200(ii,:) = min(tmp, 200000);
end
%%
figure(15)
clf
colororder({'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'});
subplot(2,2,1);
plot(tm, fact*cumsum(cumDiv_95_100,2),'linewidth',2)
title('Cumulative Diversions for 95^{th} 100 CAP')
ylabel('Water volume [km^3]')
grid on
hlegend = legend(scen.RiverName, 'Location','northwest','Orientation','horizontal');
hlegend.NumColumns=5;
hlegend.EdgeColor = 'none';
hlegend.Color = 'none';
subplot(2,2,2);
plot(tm, fact*cumsum(cumDiv_95_200,2),'linewidth',2)
title('Cumulative Diversions for 95^{th} 200 CAP')
grid on
subplot(2,2,3);
plot(tm, fact*cumsum(cumDiv_90_100,2),'linewidth',2)
title('Cumulative Diversions for 90^{th} 100 CAP')
xlabel('Time')
ylabel('Water volume [km^3]')
grid on
subplot(2,2,4);
plot(tm, fact*cumsum(cumDiv_90_200,2),'linewidth',2)
title('Cumulative Diversions for 90^{th} 200 CAP')
xlabel('Time')
grid on
% print -dpng -r600 CumulativeDiversions
%%
figure(21);
ii = 1:10;
clf
subplot(2,1,1);
plot(tm,cumDiv_90_200(ii,:)/1000,'linewidth',2)
ylabel('TAF')
if length(ii) == 1;title(scen(ii).RiverName);end
grid on
subplot(2,1,2);
plot(tm,cumsum(cumDiv_90_200(ii,:)/1000,2),'linewidth',2)
ylabel('TAF')
xlabel('Time')
grid on
%%
for ii = 1:length(c2vsim_mesh)
    a(ii,1) = polyarea(c2vsim_mesh(ii,1).X(1:end-1), c2vsim_mesh(ii,1).Y(1:end-1));
end
histogram(a*0.000247105, 30,'Normalization','probability')
%% Print the Diversions Time Series Input File
dts_fileprefix = 'DTS_DEC20_Scen';
ScenID = 1;
for prc = [95 90]
    for cap = [100 200]
        fid = fopen(fullfile('..','OptimResults', 'inputFiles',...
            [dts_fileprefix num2str(ScenID) '_' ...
            num2str(prc) '_' num2str(cap) '.dat']),'w');
        fprintf(fid, '%d %d\n', [length(scen) size(SWhyd.data,1)]);
        for ii = 1:length(scen)
            fprintf(fid,'%d ', scen(ii,1).NodeInSim);
            
            str_hyd = SWhyd.data(:,scen(ii).NodeInSim);
            pp = prctile(str_hyd, prc);
            tmp = zeros(size(SWhyd.data,1), 1);
            tmp(str_hyd > pp) = str_hyd(str_hyd > pp) - pp;
            tmp = min(tmp, cap*1000);
            tmp = tmp/1000;
            for j = 1:length(tmp)
                if tmp(j) == 0 || tmp(j) == 100  || tmp(j) == 200
                    fprintf(fid, '%.0f ', tmp(j));
                else
                    fprintf(fid, '%.3f ', tmp(j));
                end
            end
            fprintf(fid, '\n');
        end
        
        fclose(fid);
    end
end
%% Print the Diversions Element Time
fid = fopen(fullfile('..','OptimResults', 'inputFiles','divElem_DEC20_scen1.dat'),'w');
fprintf(fid,'%d\n', length(scen));
All_elem = [];
for ii = 1:length(scen)
    All_elem = [All_elem;scen(ii,1).IE];
    fprintf(fid, '%d ', [scen(ii,1).NodeInSim; length(scen(ii,1).IE); scen(ii,1).IE]);
    fprintf(fid, '\n');
end

fclose(fid);
%% Prepare the economic function
% We use the values of pred_price_county_fe column
P_LAND_TBL = readtable(fullfile('..','Rwrkspc','cv_land_price_022620.xlsx'));
% However the prices for some of the elements are missing therefore we
% interpolate the missing values
c2vsim_elem = shaperead(fullfile('..','gis_data','C2Vsim_mesh'));
p_land = nan(length(c2vsim_elem),1);
bc = zeros(length(c2vsim_elem),2);
elem_area = nan(length(c2vsim_elem),1);
for ii = 1:length(c2vsim_elem)
    bc(ii,:) = [mean(c2vsim_elem(ii,1).X(1:end-2)) mean(c2vsim_elem(ii,1).Y(1:end-2))];
    elem_area(ii,1) = polyarea(c2vsim_elem(ii,1).X(1:end-1), c2vsim_elem(ii,1).Y(1:end-1))/1e+6;
end
F = scatteredInterpolant(bc(P_LAND_TBL.ie,1), bc(P_LAND_TBL.ie,2), ...
    P_LAND_TBL.pred_price_county_fe, 'linear','nearest');

p_land(P_LAND_TBL.ie,1) = P_LAND_TBL.pred_price_county_fe;
missing_id = find(isnan(p_land));
p_land(missing_id,1) = F(bc(missing_id,1), bc(missing_id,2));
p_land(isnan(p_land),1) = 4*max(p_land);
%% Calculate the distance and the lift
x_lift = nan(length(c2vsim_elem),1);
x_distance = nan(length(c2vsim_elem),1);

for ii = 1:size(DATA,1)
    x_lift(DATA{ii,5},1) = min(0,DATA{ii,3} - DATA{ii,6});
    x_distance(DATA{ii,5},1) = DATA{ii,9}*0.621371;
end
x_lift(isnan(x_lift),1) = -10000;
x_distance(isnan(x_distance),1) = 1000;
%% Print the cost to file
fid = fopen(fullfile('..','OptimResults', 'inputFiles','ElemCost_DEC20_scen1.dat'),'w');
fprintf(fid, '%d\n', length(p_land));
fprintf(fid, '%d %.2f  %.5f  %.5f  %.5f\n',[(1:length(p_land))' elem_area p_land x_distance x_lift]');
fclose(fid);
%% Discount factors
DiscountFactors = readtable('DiscountFactors.xlsx');
fid = fopen(fullfile('..','OptimResults', 'inputFiles','DiscountFactors3.dat'),'w');
fprintf(fid,'%0.5f\n',DiscountFactors.Var2(3:end));
fclose(fid);
fid = fopen(fullfile('..','OptimResults', 'inputFiles','DiscountFactors7.dat'),'w');
fprintf(fid,'%0.5f\n',DiscountFactors.Var3(3:end));
fclose(fid);


