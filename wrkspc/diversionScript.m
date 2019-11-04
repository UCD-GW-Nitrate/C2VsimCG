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
shapewrite(c2vsim_divPoints, ['..' filesep 'gis_data' filesep 'C2Vsim_divPoints1']); 
shapewrite(c2vsim_divPolys, ['..' filesep 'gis_data' filesep 'C2Vsim_divPolys1']); 