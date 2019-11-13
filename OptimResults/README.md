In this folder you can find the optimization results.

Maximize water table rise and minimize application area
-------
This is the first optimization run that was used mainly for setting up and debbuging the tools.

We choose three diversion nodes 
1. river node 1: Kern river
2. river node 421: Kaweah river
3. river node 23: Kings river

For each diversion node we selected  a number of elements that can receive water using the existing infrastructure. (See the map)

The total number of elements that can receive diverted water is 294.

The first objective is to minimize the area.
The secod objective is to maximize the water level rise.
For the second objective we sum the (scenario - base) x weight head values for all C2Vsim nodes and all time periods. The simulation period is 1965 - 2009 with monthly step.
To encourage spatialy distributed water level rise (WLR) each meter of WLR was weighted a different value starting from 1 and reduced by 0.9.
  