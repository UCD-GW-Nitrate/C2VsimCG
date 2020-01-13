### Steps to calculate cost per elements

1.  Assign the NASS Land Use codes per mesh element.  <br />
We use an Earth Engine [Script](https://code.earthengine.google.com/4b38ae8c3362fd6da3c2eb6104648f1e) to select a year for the NASS data set. I'm currently using the year 2017. The available data span from 1997 - 2018. The script exports a GeoJSON file with a histogram object per element, which containts the distribution of NASS classes per element.

2. Data that are needed for the cost function <br />
- [land use prices](https://github.com/UCD-GW-Nitrate/C2VsimCG/blob/master/Rwrkspc/ie_c2vsim_landuse_saleprice.csv) We use the **mean_ie_lu_price** field.
- NASS - 




- Which NASS year to choose?
- What are te units in mean_ie_lu_price ($/acre)?

