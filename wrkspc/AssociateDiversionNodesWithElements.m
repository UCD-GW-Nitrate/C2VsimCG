%% Model path
c2vsim_path = ['..' filesep 'c2vsim_cg_1921ic_r374_rev' filesep 'C2VSim_CG_1921IC_R374_rev' filesep];
%% Read diversion data
[DivSpec, BypassSpec] = ReadC2VsimDiversionSpecFile([c2vsim_path 'Simulation' filesep 'CVdivspec.dat']);
river_nodes = shaperead(['..' filesep 'gis_data' filesep 'C2Vsim_riverNodes']);
river_nodes_4326 = shaperead(['..' filesep 'gis_data' filesep 'C2Vsim_riverNodes_4326']);

% make a list of canals that are not associated with a river node.
% The import water from outside the domain
ImportCanals = ListImportCanals;
%% Find the unique River IDs
IRDV_unique  = unique([DivSpec.IRDV]');
IRDV_unique(IRDV_unique == 0,:) = []; 
%% connect the unique diversion nodes with the receiving elements
clear diversionRiverNodes
diversionRiverNodes(length(IRDV_unique),1).RootName = [];
diversionRiverNodes(length(IRDV_unique),1).Names = [];
diversionRiverNodes(length(IRDV_unique),1).DivIds = [];
diversionRiverNodes(length(IRDV_unique),1).IRDV = [];
diversionRiverNodes(length(IRDV_unique),1).Coord = [];
diversionRiverNodes(length(IRDV_unique),1).Coord_4326 = [];
diversionRiverNodes(length(IRDV_unique),1).IE = [];

for ii = 1:length(IRDV_unique)
    diversionRiverNodes(ii,1).IRDV = IRDV_unique(ii);
    % find the coordinate
    ind = find([river_nodes.IRV]' == IRDV_unique(ii));
    diversionRiverNodes(ii,1).Coord = [river_nodes(ind,1).X river_nodes(ind,1).Y];
    diversionRiverNodes(ii,1).Coord_4326 = [river_nodes_4326(ind,1).X river_nodes_4326(ind,1).Y];
    inds = find([DivSpec.IRDV]' == IRDV_unique(ii));
    tmp_name = DivSpec(inds(1),1).DESCR;
    diversionRiverNodes(ii,1).Names = cell(length(inds),1);
    for j = 1:length(inds)
        diversionRiverNodes(ii,1).IE = [diversionRiverNodes(ii,1).IE; DivSpec(inds(j),1).IERELS];
        diversionRiverNodes(ii,1).Names{j,1} = DivSpec(inds(j),1).DESCR;
        diversionRiverNodes(ii,1).DivIds = [diversionRiverNodes(ii,1).DivIds;DivSpec(inds(j),1).ID];
        if j == 1
            continue;
        end
        for k = 1:length(tmp_name)
            
            if ~strcmp(tmp_name(1:k), DivSpec(inds(j),1).DESCR(1:k))
                tmp_name = DivSpec(inds(j),1).DESCR(1:k-1);
                break;
            end
        end
    end
    diversionRiverNodes(ii,1).RootName = tmp_name;
    diversionRiverNodes(ii,1).IE = unique(diversionRiverNodes(ii,1).IE);
    diversionRiverNodes(ii,1).IE(diversionRiverNodes(ii,1).IE == 0,:) = [];
end
%% Append the diversions that spread imported water
cnt_zero_divs = 0;
cnt = length(diversionRiverNodes)+1;
for ii = 1:length(ImportCanals)
    diversionRiverNodes(cnt,1).RootName = ImportCanals(ii,1).Name;
    diversionRiverNodes(cnt,1).Names = cell(0);
    diversionRiverNodes(cnt,1).IRDV = 0;
    diversionRiverNodes(cnt,1).Coord_4326 = fliplr(ImportCanals(ii,1).coords);
    n_name = length(ImportCanals(ii,1).Name);
    for jj = 1:length(DivSpec)
        if length(DivSpec(jj,1).DESCR) < n_name
            continue
        end
       if strcmp(ImportCanals(ii,1).Name, DivSpec(jj,1).DESCR(1:n_name))
           diversionRiverNodes(cnt,1).Names = [diversionRiverNodes(cnt,1).Names; DivSpec(jj,1).DESCR];
           diversionRiverNodes(cnt,1).DivIds = [diversionRiverNodes(cnt,1).DivIds; DivSpec(jj,1).ID];
           diversionRiverNodes(cnt,1).IE = [diversionRiverNodes(cnt,1).IE; DivSpec(jj,1).IERELS];
       end
    end
    diversionRiverNodes(cnt,1).IE = unique(diversionRiverNodes(cnt,1).IE);
    diversionRiverNodes(cnt,1).IE(diversionRiverNodes(cnt,1).IE == 0,:) = [];
    cnt = cnt + 1;
end

