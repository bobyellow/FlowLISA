# FlowLISA
The **univariate** local Moran's I for spatial flow data, which was proposed along with the bivariate version, namely **BiFlowLISA** (Tao and Thill, 2020). 

As detailed in Equation (3) (detailed in Tao et al. 2023), $FI_{i,j}$ is the local Moran’s I statistic, or the spatial autocorrelation measure of flow between origin $i$ and destination $j$. $f_{i,j}$ represents the value (or volume) of flow between regions $i$ and $j$. $n$ is the total number of flows in the study area. $\bar{f}$ is the average value of all flows. $w_{ij,uv}$ is the spatial flow weight between $f_{i,j}$ and $f_{u,v}$.

Equation (3):

$$
FI_{(i,j)}
\=\
\frac{
  n\\bigl(f_{(i,j)} - \bar{f}\bigr)\
  \displaystyle\sum_{(u,v)\neq(i,j)} w_{ij,uv}\\bigl(f_{(u,v)} - \bar{f}\bigr)
}{
  \displaystyle\sum_{i,j}^n \\bigl(f_{(i,j)} - \bar{f}\bigr)^2
}
$$

Run the codes of FlowLISA:
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/bobyellow/FlowLISA/blob/main/FlowLISA_main.ipynb)

To cite:

Tao, R., & Thill, J. C. (2020). BiFlowLISA: Measuring spatial association for bivariate flow data. Computers, Environment and Urban Systems, 83, 101519.

Tao, R., Chen, Y., & Thill, J. C. (2023). A space-time flow LISA approach for panel flow data. Computers, Environment and Urban Systems, 106, 102042.
