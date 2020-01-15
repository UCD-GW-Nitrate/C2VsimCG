### Steps to calculate cost per element

1.  Assign the NASS Land Use codes per mesh element.  <br />
We use an Earth Engine [Script](https://code.earthengine.google.com/4b38ae8c3362fd6da3c2eb6104648f1e) (this is private and I haven't found a way to share to a non earth engine user) to select a year for the NASS data set. I'm currently using the year 2017. The available data in earth engine span from 1997 - 2018. The script exports a GeoJSON file with a histogram object per finite element, which containts the distribution of NASS classes per element.

2. Data that are needed for the cost function <br />
- [land use prices](https://github.com/UCD-GW-Nitrate/C2VsimCG/blob/master/Rwrkspc/ie_c2vsim_landuse_saleprice.csv) We use the **mean_ie_lu_price** field.
- [NASS - Model association](https://github.com/UCD-GW-Nitrate/C2VsimCG/blob/master/Rwrkspc/ca_cdl_xwalk.xlsx) This makes a correspondance between the NASS codes and the class used by the model. 
- Depth to water table. This is calculated using the average depth of a selected time period.

3. Cost function logic <br />
For each element we have a list of NASS codes along with coverage percentages.<br />
Then we loop through each NASS crop code. <br />
We first identify the model class (davis_class) using the association file.<br />
Next we get a list of all elements that have a price for that davis class. We distinguish 3 cases: 
<br /> 1. There is a price listed for that crop in the current element.
<br /> 2. There is no price in the current element for the crop type but there are more that one elements with prices within a short distance from the element. Then we select the element that the simulated water table depth is closest to the element we are currently calculate the cost function.
<br /> 3. There is only one or no element in short distance with price. Then we use the element which is closer to the current element. In some cases that distance was quite significant.
We computed the area weighted average for each element.

## Results
Preliminary results can be found at https://bit.ly/35OxMhq

#### Clarifications
- Which NASS year to choose (2017)?
- What are te units in mean_ie_lu_price ($/acre)?
- Which period to use for depth to water table (average of the last decade)?
- What is a short distance (20 km)?

