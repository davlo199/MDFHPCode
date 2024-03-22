The data files can be found publicly from [USGS quake search](https://earthquake.usgs.gov/earthquakes/search/). 
We are only concerned with the columns 1 (A) called "time" and 5 (E) called "mag" containing the occurrence date of the events and magnitude specifically. 

Data cleaning i.e. turning date time into a numeric value from the origin time and selecting the correct magnitude and time window is done at the start of every script that reads in the data.
