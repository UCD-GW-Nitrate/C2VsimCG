%% For the objective function we would need the timeseries
% As first attempt we will divert water from 3 stream locations only.
% run the snipet in the diversionScrpit that generates the variable D
%
% The 3 diversion points are the ones with ids 1, 3, 4
% Calculate the time series
load('c2vsimBasicScenario.mat', 'SWhyd')

divID = [];
divTimeSeries = [];
for ii = [1 3 4]
    [c, d] = max(sum(SWhyd.data(:,D(ii,1).Nodes)));
    Qt = SWhyd.data(:,D(ii,1).Nodes(d));
    [f,x] = ecdf(Qt);
    q = interp1(f,x,0.95);
    Qt(Qt < q) = 0;
    divID = [divID; D(ii,1).Nodes(d)];
    divTimeSeries = [divTimeSeries Qt];
end
%% write into file
fid = fopen(['..' filesep 'RunC2Vsim' filesep 'divTimeSeries.dat'],'w');
fprintf(fid, '%d %d\n' , [length(divID) size(divTimeSeries,1)]);
for ii = 1:length(divID)
    fprintf(fid, '%d', divID(ii));
    fprintf(fid, ' %f', divTimeSeries(:,ii)/1000);
    fprintf(fid, '\n');
end
fclose(fid);
%% 
load('c2vsimDiversionData.mat', 'divElems');
%% write diversions-elements to file
id = find(isnan([divElems.DivNode]'));
divElems(id,:) = [];
divID = unique([divElems.DivNode]');

fid = fopen(['..' filesep 'RunC2Vsim' filesep 'divElem.dat'],'w');
fprintf(fid, '%d\n', length(divID));
for ii = 1:length(divID)
    id = find([divElems.DivNode]' == divID(ii));
    fprintf(fid, '%d %d', [divID(ii) length(id)]);
    fprintf(fid, ' %d', [divElems(id,1).IE]);
    fprintf(fid, '\n');
end

fclose(fid);
%% write element info file
load('c2vsimGeometryData.mat', 'c2vsim_mesh')
for ii = 1:size(c2vsim_mesh,1)
    area(ii,1) = polyarea(c2vsim_mesh(ii,1).X(1:end-1), c2vsim_mesh(ii,1).Y(1:end-1))/(1000*1000);
end
fid = fopen(['..' filesep 'RunC2Vsim' filesep 'ElemInfo.dat'],'w');
fprintf(fid,'%d\n', size(c2vsim_mesh,1));
fprintf(fid, '%d %f\n', [[c2vsim_mesh.IE]' area]');
fclose(fid);
