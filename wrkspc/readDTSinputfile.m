function DTS = readDTSinputfile(filename)
% readDTSinputfile Reads the input file of the optimization
% that describe % the diversion time series. 
% e.g the node ids and the diversion time series

fid = fopen(filename);
cc = textscan(fid, '%f %f',1);
for ii = 1:cc{1,1}
    ccc = textscan(fid, '%f',cc{1,2} + 1);
    DTS(ii,1).Id = ccc{1,1}(1);
    DTS(ii,1).Data = ccc{1,1}(2:end);
end

fclose(fid);