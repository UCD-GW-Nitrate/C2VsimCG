%% Extrach streamhydrographs for all nodes
% Run this after modifing the CVprint.dat file to print all streams
SWhyd = readC2VsimSWhyd('../RunC2Vsim/Results/CVSWhyd.out');
%% save
save('c2vsimBasicScenario', 'SWhyd');