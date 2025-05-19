# FlowLISA
The **univariate** local Moran's I for spatial flow data, which was proposed along with the bivariate version, namely **BiFlowLISA**, in the paper:
Tao, R., & Thill, J. C. (2020). BiFlowLISA: Measuring spatial association for bivariate flow data. Computers, Environment and Urban Systems, 83, 101519.

Tao and Thill (2020) extended the local Moran’s I to spatial flow data.

As detailed in Equation (3), $FI_{i,j}$ is the local Moran’s I statistic, or the spatial autocorrelation measure of flow between origin $i$ and destination $j$; $f_{i,j}$ represents the value (or volume) of flow between regions $i$ and $j$; $n$ is the total number of flows in the study area; $f$ is the average value of all flows; $w_{ij,uv}$ is the spatial flow weight between $f_{i,j}$ and $f_{u,v}$.

Equation (3):
$$
FI_{i,j} = \frac{n f_{i,j} - \sum_{u,v \neq i,j} w_{ij,uv}f_{u,v} - f}{\sum_{i,j} n (f_{i,j} - f)^2}
$$

Run the codes of FlowLISA:
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bobyellow/FlowLISA/blob/main/FlowLISA_main.ipynb)
