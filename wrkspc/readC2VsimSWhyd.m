function SWhyd = readC2VsimSWhyd(filename)
fid = fopen(filename, 'r');
%read 5 header files
for i = 1:5; temp = fgetl(fid);end

%read stream nodes
temp = fgetl(fid);
c = strsplit(temp, {' ','\t'});
SWhyd.Nodes = cell2mat(cellfun(@str2double,c,'UniformOutput',false))';


cnt_per = 0;
while 1
    temp = fgetl(fid);
    if temp == -1
        break;
    end
    if isempty(temp)
        break;
    end
    C = strsplit(temp, {' ','\t'});
    c = textscan(C{1,1},'%f/%f/%f_%s');
    cnt_per = cnt_per + 1;
    SWhyd.time{cnt_per,1} = [num2str(c{1,1}) '/' num2str(c{1,2}) '/' num2str(c{1,3})];
    tmp = cell2mat(cellfun(@str2double,C,'UniformOutput',false))';
    tmp(isnan(tmp),:) = [];
    SWhyd.data(cnt_per, :) = tmp;
    %display(SWhyd.time{cnt_per,1})
    %if exist ('OCTAVE_VERSION', 'builtin') > 0
    %    fflush(stdout);
    %end
    %for j = 1:Nswhyd
    %    SWhyd.data(cnt_per, j) = str2double(C{1,1+j});
    %end
end
fclose(fid);