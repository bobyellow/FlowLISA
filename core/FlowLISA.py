import numpy as np
from core.getFlowNeighbors import getFlowNeighborsContiguity  # Obtain neighbors of flow data
from core.spatstats import (
    calculateMoranI, 
    calculateGetisG, 
    calculateGearyC, 
    calculateMultiGearyC,
    
)


# Define the univariate flow LISA of which users can select the foundational spatial stat

def execFLOWLISA(AREAS1, AREAS2, FlowValue, Spatstat, NeiLvl):
    """
    Execute FlowLISA to analyze spatial autocorrelation in univariate flow data.

    Parameters:
    1. AREAS1: Origin areas (list of polygons)
    2. AREAS2: Destination areas (list of polygons)
    3. FlowValue: Dictionary of (O, D) flow values
    4. Spatstat:
        1 -> Local Moran's I
        2 -> Local Getis-Ord G
        3 -> Local Geary's C
    5. NeiLvl: Neighborhood level for flow connections

    Returns:
    - A formatted output string containing results
    """
    
    areas1, areas2 = AREAS1, AREAS2  # Assign origin and destination areas
    y = FlowValue  # Dictionary containing flow values
    yKeys, yValues = list(y.keys()), list(y.values())

    # Initialize output dictionary with original values
    yOutput = {k: [v] if Spatstat != 5 else v for k, v in y.items()}

    # Obtain flow neighbors based on O&D contiguity (default: Rook's case)
    Wflow = getFlowNeighborsContiguity(areas1, areas2, y, NeiLvl)

    # Compute global statistics for local calculations
    dataSum, dataMean, dataStd = np.sum(yValues), np.mean(yValues), np.std(yValues)
    GMoranI = 0  # Global Moran's I initialization

    # Calculate spatial statistics for each flow
    for s in yKeys:
        neighbors = Wflow.get(s, [])
        
        if Spatstat == 1 and neighbors:
            MoranI = calculateMoranI(s, neighbors, dataMean, dataStd, y, len(yKeys))
            #print(MoranI)
            yOutput[s].extend([MoranI, 0])  # Moran's I and initial p-value
            GMoranI += MoranI
        elif Spatstat == 2:
            GetisG = calculateGetisG(neighbors, dataMean, dataStd, y, len(yKeys)) if neighbors else 0
            #print(GetisG)
            yOutput[s].extend([GetisG, 0])
        elif Spatstat == 3:
            GearyC = calculateGearyC(s, neighbors, y) if neighbors else 0
            yOutput[s].extend([GearyC, 0])


    # Monte-Carlo significance testing (1000 simulations)
    GMoranI_sim = np.zeros(1000)

    for i in range(1000):
        # Ensure randomized keys are tuples
        shuffled_keys = [tuple(key) for key in np.random.permutation(yKeys)]  # Convert to list of tuples

        # Create randomized dictionary
        yRandom = {rk: y[tuple(ok)] for rk, ok in zip(shuffled_keys, yKeys)}
        
        for pk in yRandom:
            neighbors = Wflow.get(pk, [])
            if Spatstat == 1 and neighbors:
                MoranI = calculateMoranI(pk, neighbors, dataMean, dataStd, yRandom, len(yKeys))
                GMoranI_sim[i] += MoranI
                if abs(MoranI) > abs(yOutput[pk][1]):
                    yOutput[pk][2] += 1
            elif Spatstat == 2 and neighbors:
                GetisG = calculateGetisG(neighbors, dataMean, dataStd, yRandom, len(yKeys))
                if GetisG > yOutput[pk][1]:
                    yOutput[pk][2] += 1
            elif Spatstat == 3 and neighbors:
                GearyC = calculateGearyC(pk, neighbors, yRandom) 
                if GearyC > yOutput[pk][1]:
                    yOutput[pk][2] += 1


    # Analyze Monte-Carlo simulation results
    GMoranI_sim.sort()
    significance_msg = (
        f"Global Moran's I value is: {GMoranI:.4f}. It is "
        f"{'positive' if GMoranI >= 0 else 'negative'}, "
        f"{'significantly' if GMoranI >= GMoranI_sim[950] or GMoranI <= GMoranI_sim[49] else 'insignificantly'} at 0.01 level."
    )

    # Construct output strings based on Spatstat type
    output_str_list = []
    if Spatstat == 1:
        output_str_list.append(significance_msg)
        output_str_list.append('O, D, V, MoranI, p-value, pattern')
        for k, v in yOutput.items():
            v[2] /= 1000.0
            significance = "NS"
            if v[1] != 0 or v[2] != 0:
                if v[2] <= 0.05:
                    significance = "HH" if v[0] > dataMean else "LL" if v[1] > 0 else "HL" if v[0] > dataMean else "LH"
            output_str_list.append(f"{k[0]}, {k[1]}, {', '.join(map(str, v))}, {significance}")
    
    if Spatstat == 2:
        output_str_list.append('O, D, V, GetisG, p-value, pattern')
        for k, v in yOutput.items():
            v[2] /= 1000.0
            significance = "NS"
            if v[1] != 0 or v[2] != 0:
                if v[2] <= 0.05:
                    significance = "hot" 
                elif v[2] >= 0.95:
                    significance = "cold" 
            output_str_list.append(f"{k[0]}, {k[1]}, {', '.join(map(str, v))}, {significance}")

    if Spatstat == 3:
        output_str_list.append('O, D, V, GearyC, p-value, pattern')
        for k, v in yOutput.items():
            v[2] /= 1000.0
            significance = "NS"
            if v[1] != 0 or v[2] != 0:
                if v[2] <= 0.05:
                    significance = "dissimilar" 
                elif v[2] >= 0.95:
                    significance = "similar" 
            output_str_list.append(f"{k[0]}, {k[1]}, {', '.join(map(str, v))}, {significance}")

    return '\n'.join(output_str_list)





