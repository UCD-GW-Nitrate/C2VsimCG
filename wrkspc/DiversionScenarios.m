function scen = DiversionScenarios(ScenarioName)
% RiverName is the name of the river
% NodeInSim is the river ID where the water is extracted from during the
%               optimization
% DivIds are the indices in the diversionRiverNodes structure that is
% return from the AssociateDiversionNodesWithElements script
% ImportNode
% rivID is the id of the river in the river line shapefile C2Vsim_rivers_3310.shp

switch ScenarioName
    case 'Scen1'
        scen(1,1).RiverName = 'Calaveras';
        scen(1,1).NodeInSim = 166;
        scen(1,1).DivIds = 24;
        scen(1,1).ImportNode = 161;
        scen(1,1).rivID = 25;
        
        scen(2,1).RiverName = 'Stanislaus';
        scen(2,1).NodeInSim = 153;
        scen(2,1).DivIds = [23; 86];
        scen(2,1).ImportNode = 146;
        scen(2,1).rivID = 23;
        
        scen(3,1).RiverName = 'Tuolumne';
        scen(3,1).NodeInSim = 141;
        scen(3,1).DivIds = [20; 21];
        scen(3,1).ImportNode = 135;
        scen(3,1).rivID = 21;
        
        scen(4,1).RiverName = 'Merced';
        scen(4,1).NodeInSim = 120;
        scen(4,1).DivIds = [17; 18; 83];
        scen(4,1).ImportNode = 116;
        scen(4,1).rivID = 17;
        
        scen(5,1).RiverName = 'Chowchilla';
        scen(5,1).NodeInSim = 83;
        scen(5,1).DivIds = [15];
        scen(5,1).ImportNode = 80;
        scen(5,1).rivID = 11;
        
        scen(6,1).RiverName = 'Fresno';
        scen(6,1).NodeInSim = 72;
        scen(6,1).DivIds = [14];
        scen(6,1).ImportNode = 69;
        scen(6,1).rivID = 9;
        
        scen(7,1).RiverName = 'San Joaquin A';
        scen(7,1).NodeInSim = 61;
        scen(7,1).DivIds = [13];
        scen(7,1).ImportNode = 54;
        scen(7,1).rivID = 7;
        
        scen(8,1).RiverName = 'San Joaquin B';
        scen(8,1).NodeInSim = 156;
        scen(8,1).DivIds = [16;19;22];
        scen(8,1).ImportNode = 134;
        scen(8,1).rivID = [16; 18; 20; 22; 24];
        
        scen(9,1).RiverName = 'Kings';
        scen(9,1).NodeInSim = 29;
        scen(9,1).DivIds = [7; 8; 9; 10; 11; 12];
        scen(9,1).ImportNode = 23;
        scen(9,1).rivID = 2;
        
        scen(10,1).RiverName = 'Kaweah';
        scen(10,1).NodeInSim = 427;
        scen(10,1).DivIds = [72; 73; 74; 75];
        scen(10,1).ImportNode = 420;
        scen(10,1).rivID = 5;
        
        scen(11,1).RiverName = 'Tule';
        scen(11,1).NodeInSim = 19;
        scen(11,1).DivIds = [5; 6];
        scen(11,1).ImportNode = 10;
        scen(11,1).rivID = 6;
        
        scen(12,1).RiverName = 'Kern';
        scen(12,1).NodeInSim = 5;
        scen(12,1).DivIds = [1; 2; 3; 4];
        scen(12,1).ImportNode = 1;
        scen(12,1).rivID = 1;
        
    case 'Scen2'
        scen = DiversionScenarios('Scen1');
        scen(13,1).RiverName = 'Friant-Kern Canal';
        scen(13,1).NodeInSim = 54;
        scen(13,1).DivIds = 77;
        scen(13,1).ImportNode = 54;
        scen(13,1).rivID = [];
    case 'Scen3'
        scen = DiversionScenarios('Scen2');
        scen(14,1).RiverName = 'Delta';
        scen(14,1).NodeInSim = 418;
        scen(14,1).DivIds = [71; 79; 80];
        scen(14,1).ImportNode = 418;
        scen(14,1).rivID = [];
    otherwise
        scen = [];
end

