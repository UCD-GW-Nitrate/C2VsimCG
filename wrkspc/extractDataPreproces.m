
%% C2Vsim path
c2vsim_path = ['..' filesep 'c2vsim_cg_1921ic_r374_rev' filesep 'C2VSim_CG_1921IC_R374_rev' filesep];
%% Read Mesh nodes
fid = fopen([c2vsim_path 'Preprocessor' filesep 'CVnode.dat'],'r');
temp = textscan(fid, '%f %f %f', 1393, 'HeaderLines',80);
fclose(fid);
XY = [temp{1,2} temp{1,3}];
ND_ID = temp{1,1};
%% Read Mesh Elements
fid = fopen([c2vsim_path 'Preprocessor' filesep 'CVelement.dat'],'r');
temp = textscan(fid, '%f %f %f %f %f', 1392, 'HeaderLines',93);
fclose(fid);
MSH = [temp{1,2} temp{1,3} temp{1,4} temp{1,5}];
EL_ID = temp{1,1};
%% Read Element Characteristics file
fid = fopen([c2vsim_path 'Preprocessor' filesep 'CVcharac.dat'],'r');
CHARAC = textscan(fid, '%f %f %f %f %f %f %f', 1392, 'HeaderLines',78);
fclose(fid);
%% Read Stratigraphy file
fid = fopen([c2vsim_path 'Preprocessor' filesep 'CVstrat.dat'],'r');
strat = textscan(fid, '%f %f %f %f %f %f %f %f %f', 1393, 'HeaderLines',92);
fclose(fid);
%% Define lake elements from CVlake.dat
lake_elements = [1352;1353;1363;1364;1109;1110;1111;1136;1137;1138];
%% Write Element shapefile
clear S
for i = 1:size(MSH,1)
    if MSH(i,4) == 0
        nd = [1 2 3 1];
    else
        nd = [1 2 3 4 1];
    end
    xx = XY(MSH(i,nd),1);
    yy = XY(MSH(i,nd),2);
    
    S(i,1).Geometry = 'Polygon';
    S(i,1).BoundingBox = [min(xx) min(yy); ...
                          max(xx) max(yy)];
    S(i,1).X = [xx' nan];
    S(i,1).Y = [yy' nan];
    S(i,1).id = 1;
    S(i,1).IE = EL_ID(i);
    S(i,1).IRNE = CHARAC{1,2}(i);
    S(i,1).FRNE = CHARAC{1,3}(i);
    S(i,1).ISTE = CHARAC{1,4}(i);
    S(i,1).IRGE = CHARAC{1,5}(i);
    S(i,1).ISGE = CHARAC{1,6}(i);
    S(i,1).ISOILE = CHARAC{1,7}(i);
    S(i,1).is_lake = 0;
    if find(lake_elements == i)
        S(i,1).is_lake = 1;
    end
end
c2vsim_mesh = S;
%% Write Mesh nodes
clear S
for i = 1:size(XY,1)
    xx = XY(i,1);
    yy = XY(i,2);
    
    S(i,1).Geometry = 'Point';
    S(i,1).X = xx;
    S(i,1).Y = yy;
    S(i,1).BoundingBox = [min(xx) min(yy); ...
                          max(xx) max(yy)];
    S(i,1).ID = ND_ID(i);
    S(i,1).ELV = strat{1,2}(i);
    S(i,1).W1 = strat{1,3}(i);
    S(i,1).W2 = strat{1,4}(i);
    S(i,1).W3 = strat{1,5}(i);
    S(i,1).W4 = strat{1,6}(i);
    S(i,1).W5 = strat{1,7}(i);
    S(i,1).W6 = strat{1,8}(i);
end
c2vsim_nodes = S;

%% Extract Rivers
fid = fopen([c2vsim_path 'Preprocessor' filesep 'CVrivers.dat'],'r');
% read NRH, NR, NRTB
temp = textscan(fid, '%f / %s', 3, 'HeaderLines',64);
NRH = temp{1,1}(1); %NUmber of stream reaches
NR = temp{1,1}(2); %Number of stream nodes
NRTB = temp{1,1}(3); %Number of data points in tables per stream node

% read comments
for ii = 1:22; tline = fgetl(fid);end

for i = 1:NRH
    i
    % for each stream read
    % a comment line
    while 1
        tline = fgetl(fid);
        if strcmp('C--------------------------------------------',tline)
            break;
        end
    end
    % the name of the stream
    temp = textscan(fid, 'C     REACH  %f  -  %[^\n]', 1);
    %temp = textscan(fid, 'C     REACH  %f  -  %s', 1);
    RIVER(i,1).name = temp{1,2}{1,1}(1:end-1);
    % next 4 comment lines
    for ii = 1:5; tline = fgetl(fid);end
    temp = textscan(fid, '%f %f %f %f', 1);
    ID = temp{1,1}(1);
    IBUR = temp{1,2}(1);
    IBDR = temp{1,3}(1);
    IDWN = temp{1,4}(1);
    RIVER(i,1).ID = ID;
    RIVER(i,1).IBUR = IBUR;
    RIVER(i,1).IBDR = IBDR;
    RIVER(i,1).IDWN = IDWN;
    % 5 more comment lines
    for ii = 1:6; tline = fgetl(fid);end
    temp = textscan(fid, '%f %f %f', IBDR-IBUR+1);
    RIVER(i,1).IRV = temp{1,1};
    RIVER(i,1).IGW = temp{1,2};
    RIVER(i,1).IRGST = temp{1,3};
    tline = fgetl(fid);
end

% read comments and other lines
for ii = 1:35; tline = fgetl(fid);end

% read stream rate tables
for i = 1:NR
    temp = textscan(fid, '%f %f %f %f', 1);
    StreamTable(i,1).ID = temp{1,1};
    StreamTable(i,1).BOTR = temp{1,2};
    StreamTable(i,1).HRTB = temp{1,3};
    StreamTable(i,1).QRTB = temp{1,3};
    temp = textscan(fid, '%f %f', 4);
    StreamTable(i,1).HRTB = [StreamTable(i,1).HRTB; temp{1,1}];
    StreamTable(i,1).QRTB = [StreamTable(i,1).QRTB; temp{1,2}];
end
fclose(fid);
%% Write Rivers shapefile
clear S
for i = 1:NRH
    xx = XY(RIVER(i,1).IGW,1);
    yy = XY(RIVER(i,1).IGW,2);
    S(i,1).Geometry = 'Polyline';
    S(i,1).BoundingBox = [min(xx) min(yy); ...
                          max(xx) max(yy)];
    S(i,1).X = [xx' nan];
    S(i,1).Y = [yy' nan];
    S(i,1).ID = RIVER(i,1).ID;
    S(i,1).name = RIVER(i,1).name;
end
c2vsim_rivers = S;
%% Write River nodes
clear S
% Define inflows nodes
inflows = [205;211;220;218;225;233;243;237;248;256;263;269;283;341;349;357;390;374;400;188;182;173;161;146;135;128;116;105;93;80;69;54;23;420;10;1;24;11;421;4;4];
% make a list of nodes
id_list = [];
for i = 1:size(RIVER,1)
    for j = 1:size(RIVER(i,1).IRV,1)
        id_list = [id_list; RIVER(i,1).IRV(j) RIVER(i,1).IGW(j)];
    end
end
for i = 1:size(id_list,1)
    xx = XY(id_list(i,2),1);
    yy = XY(id_list(i,2),2);
    S(i,1).Geometry = 'Point';
    S(i,1).BoundingBox = [min(xx) min(yy); ...
                          max(xx) max(yy)];
    S(i,1).X = xx;
    S(i,1).Y = yy;
    S(i,1).IGW = id_list(i,2);
    S(i,1).IRV = id_list(i,1);
    S(i,1).Inflow = 0;
    if ~isempty(find(inflows == S(i,1).IRV))
        S(i,1).Inflow = 1;
    end
end
c2vsim_riverNodes = S;
%% Read wells
fid = fopen([c2vsim_path 'Preprocessor' filesep 'CVWells.dat'],'r');
wells = textscan(fid, '%f %f %f %f %f %f / %[^\n]', 133, 'HeaderLines',78);
fclose(fid);
%% Write wells shapefile
clear S
for i = 1:133
    xx = wells{1,2}(i);
    yy = wells{1,3}(i);
    S(i,1).Geometry = 'Point';
    S(i,1).BoundingBox = [min(xx) min(yy); ...
                          max(xx) max(yy)];
    S(i,1).X = xx;
    S(i,1).Y = yy;
    S(i,1).ID = wells{1,1}(i);
    S(i,1).name = wells{1,7}{i};
    S(i,1).RWELL = wells{1,4}(i);
    S(i,1).PERFT = wells{1,5}(i);
    S(i,1).PERFB = wells{1,6}(i);
end
c2vsim_wells = S;
%% Read small watershed boundaries
% !!!!to read this remove the lines 920-922 from the file but then paste them back again!!!!!
fid = fopen([c2vsim_path 'Simulation' filesep 'CVbound.dat'],'r');
C = textscan(fid,'%s',1, 'headerlines',518);
for i = 1:210
	C = textscan(fid,'%d %d %d %d %d %d %d',1);
	SWBND(i,1).ID = C{1,1};
	SWBND(i,1).IWBS = C{1,2};
	SWBND(i,1).AREAS = C{1,3};
	SWBND(i,1).IWBTS = C{1,4};
	SWBND(i,1).IWB = C{1,6};
	SWBND(i,1).QMAXWB = C{1,7};
	for j = 1:C{1,5}-1
		C = textscan(fid,'%d %d',1);
        SWBND(i,1).IWB = [SWBND(i,1).IWB C{1,1}];
        SWBND(i,1).QMAXWB = [SWBND(i,1).QMAXWB C{1,2}];
    end
end
fclose(fid);
%%  Create shapefile with the Small watershed segments
clear S;
for i = 1:length(SWBND)
    S(i,1).Geometry = 'PolyLine';
    S(i,1).BoundingBox = zeros(2);
    S(i,1).X = [XY(SWBND(i,1).IWB,1)' nan];
    S(i,1).Y = [XY(SWBND(i,1).IWB,2)' nan];
    S(i,1).BoundingBox = [nanmin(S(i,1).X) nanmin(S(i,1).Y); ...
                          nanmax(S(i,1).X) nanmax(S(i,1).Y)];
    S(i,1).ID = double(SWBND(i,1).ID);
    S(i,1).IWBS = double(SWBND(i,1).IWBS);
    S(i,1).AREAS = double(SWBND(i,1).AREAS);
    S(i,1).IWBTS = double(SWBND(i,1).IWBTS);
    if length(SWBND(i,1).QMAXWB) == 1
        S(i,1).QMAXWB = double(SWBND(i,1).QMAXWB(1));
    else
        S(i,1).QMAXWB = double(SWBND(i,1).QMAXWB(2));
    end
end
c2vsim_wshed = S;
%% create a shapefile with the nodes of each stream watershed
cnt = 1;
clear S;
for i = 1:length(c2vsim_wshed)
    for j = 1:length(c2vsim_wshed(i,1).X)-1
        xx = c2vsim_wshed(i,1).X(j);
        yy = c2vsim_wshed(i,1).Y(j);
        S(cnt,1).Geometry = 'Point';
        S(cnt,1).BoundingBox = [min(xx) min(yy); ...
                              max(xx) max(yy)];
        S(cnt,1).X = xx;
        S(cnt,1).Y = yy;
        S(cnt,1).WSHD_ID = c2vsim_wshed(i,1).ID;
        cnt = cnt + 1;
    end
end
c2vsim_wshedNodes = S;
%% Save data
save('c2vsimGeometryData','c2vsim_mesh', 'c2vsim_nodes', 'c2vsim_rivers',...
    'c2vsim_riverNodes', 'c2vsim_wells', 'c2vsim_wshed', 'c2vsim_wshedNodes');
%% write data to shapefiles
load('c2vsimGeometryData.mat');
shapewrite(c2vsim_mesh, ['..' filesep 'gis_data' filesep 'C2Vsim_mesh1']); 
shapewrite(c2vsim_nodes, ['..' filesep 'gis_data' filesep 'C2Vsim_nodes1']); 
shapewrite(c2vsim_rivers, ['..' filesep 'gis_data' filesep 'C2Vsim_rivers1']); 
shapewrite(c2vsim_riverNodes, ['..' filesep 'gis_data' filesep 'C2Vsim_riverNodes1']);
shapewrite(c2vsim_wells, ['..' filesep 'gis_data' filesep 'C2Vsim_wells1']);
shapewrite(c2vsim_wshed, ['..' filesep 'gis_data' filesep 'C2Vsim_wshed1']);
shapewrite(c2vsim_wshedNodes, ['..' filesep 'gis_data' filesep 'C2Vsim_wshedNodes1']);