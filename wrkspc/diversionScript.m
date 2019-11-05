%% Paths
c2vsim_path = ['..' filesep 'c2vsim_cg_1921ic_r374_rev' filesep 'C2VSim_CG_1921IC_R374_rev' filesep];
%%
[DivSpec, BypassSpec] = ReadC2VsimDiversionSpecFile([c2vsim_path 'Simulation' filesep 'CVdivspec.dat']);
%%
Diversions = ReadC2VsimDiversionsFile([c2vsim_path 'Simulation' filesep 'CVdiversions.dat']);
%% make a shape file with the diversion points that divert water from streams only
load('c2vsimGeometryData', 'c2vsim_riverNodes')
clear S
divIds = [DivSpec.IRDV]';
divIds = unique(divIds(divIds ~= 0));
for ii = 1:length(divIds)
    iriv = find([c2vsim_riverNodes.IRV]' == divIds(ii));
    xx = c2vsim_riverNodes(iriv,1).X;
    yy = c2vsim_riverNodes(iriv,1).Y;
    S(ii,1).Geometry = 'Point';
    S(ii,1).BoundingBox = [min(xx) min(yy); ...
                           max(xx) max(yy)];
    S(ii,1).X = xx;
    S(ii,1).Y = yy;
    S(ii,1).IRV = divIds(ii);
end
c2vsim_divPoints = S;
%% make a shapefiles with the polygons that receive diversions
clear S
load('c2vsimGeometryData', 'c2vsim_mesh')
polyids = [];
for ii = 1:size(DivSpec,1)
    if DivSpec(ii,1).IRDV == 0
        continue;
    end
    polyids = [polyids; DivSpec(ii,1).IERELS];
end
polyids = unique(polyids);
polyids(polyids==0,:) = [];
c2vsim_divPolys = c2vsim_mesh(polyids,1);
%%
save('c2vsimDiversionData','DivSpec', 'BypassSpec', 'Diversions', 'c2vsim_divPoints','c2vsim_divPolys')
%% Write to shapefiles
load('c2vsimDiversionData','c2vsim_divPoints','c2vsim_divPolys');
shapewrite(c2vsim_divPoints, ['..' filesep 'gis_data' filesep 'C2Vsim_divPoints']); 
shapewrite(c2vsim_divPolys, ['..' filesep 'gis_data' filesep 'C2Vsim_divPolys']); 
%%
D(1,1).Name = 'KERN RIVER';
D(1,1).Nodes = 1:9;
D(2,1).Name = 'TULE RIVER';
D(2,1).Nodes = 10:22;
D(3,1).Name = 'KAWEAH RIVER';
D(3,1).Nodes = 420:428;
D(4,1).Name = 'KINGS RIVER';
D(4,1).Nodes = 23:32;
D(5,1).Name = 'SAN JOAQUIN RIVER';
D(5,1).Nodes = 54:62;
D(6,1).Name = 'FRESNO RIVER';
D(6,1).Nodes = 69:77;
D(7,1).Name = 'CHOWCHILLA RIVER';
D(7,1).Nodes = 80:88;
D(8,1).Name = 'DEADMAN''S'' CREEK';
D(8,1).Nodes = 93:101;
D(9,1).Name = 'BEAR CREEK';
D(9,1).Nodes = 105:112;
D(10,1).Name = 'MERCED RIVER';
D(10,1).Nodes = 116:124;
D(11,1).Name = 'TUOLUMNE RIVER';
D(11,1).Nodes = 135:143;
D(12,1).Name = 'STANISLAUS RIVER';
D(12,1).Nodes = 146:154;
D(13,1).Name = 'CALAVERAS RIVER';
D(13,1).Nodes = 161:167;
%% plot cumulative stream hydrograph
ii = 3;
plot(cumsum(SWhyd.data(:,D(ii,1).Nodes)));
[c, d] = max(sum(SWhyd.data(:,D(ii,1).Nodes)));
title([D(ii,1).Name ', #: ' num2str(D(ii,1).Nodes(d))])
%% plot ECDF of the stream hydrograph
ii = 4;
[c, d] = max(sum(SWhyd.data(:,D(ii,1).Nodes)));
Qt = SWhyd.data(:,D(ii,1).Nodes(d));
[f,x] = ecdf(Qt);
q = interp1(f,x,0.95);
Qt(Qt < q) = 0;
plot(Qt)
title(num2str(sum(Qt)))
