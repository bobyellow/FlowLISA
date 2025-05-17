import core.shapefile
import pandas as pd
from collections import defaultdict
from core.FlowLISA import execFLOWLISA

# Import input flow data from .txt files
flowdf1 = pd.read_csv('input/Flow37xLL.txt', sep='\s+')
F_dt1 = dict(zip(zip(flowdf1['O'], flowdf1['D']), flowdf1['Flow']))

# The input flow data should not contain zero-value flow (OD pair)
# The input flow data are stored as dictionary format, (O,D) tuple as key, flow values as lists

# Import Origin and Destination shapefiles using core.shapefile
StationPolygon1 = core.shapefile.Reader("input/Hex37_O.shp")
StationPolygon2 = core.shapefile.Reader("input/Hex37_D.shp")

# Extract polygon shapes
shapes1 = StationPolygon1.shapes()
shapes2 = StationPolygon2.shapes()

# Prepare AREAS input for Queen's and Rook's contiguity
AREAS1 = [[shape.points] for shape in shapes1]  # Ensure proper structure for AREAS
AREAS2 = [[shape.points] for shape in shapes2]  # Ensure proper structure for AREAS

# Execute FlowLISA function
outputStr = execFLOWLISA(AREAS1, AREAS2, F_dt1, 1, 120)
"""
    Execute FlowLISA to analyze spatial autocorrelation in univariate flow data
    
    Parameters of execFLOWLISA(AREAS1, AREAS2, FlowValue, Spatstat, NeiLvl):
    1. AREAS1: Origin areas (list of polygons)
    2. AREAS2: Destination areas (list of polygons)
    3. FlowValue: Dictionary of (O, D) flow values
    4. Spatstat:
        1 -> Local Moran's I
        2 -> Local Getis-Ord G
        3 -> Local Geary's C
    5. NeiLvl: Neighborhood level for flow connections
        # Level=1: one of OD is the same and the other is neighbor;
        # Level=2: both OD are neighbors, so level ==12 means a combination of the two above
        # 18 means same D, Os are neighbors
        # 19 means same O, Ds are neighbors
        # adding 0 means including the situation of flows sharing the same O & D as flow i
        # refer to getFlowNeighbors.py for more details

    Returns:
    - A formatted output string containing results
"""
   

# Save output to text file
output_filename = 'result/FlowLISA_I_Fake37xLL_Nei120_0001.txt'
with open(output_filename, 'w') as outputFile:
    outputFile.write(outputStr)

print(f"Processing complete. Results saved to {output_filename}")
