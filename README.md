# Climate Data-Analysis
Visualization
This repository contains code and documentation for analyzing climate data. 
The following outline aims to provide a brief description of the project's goals and objectives.

## Table of Contents

- [Introduction]
- [Reading Netcdf]
- [Dataset Analysis]
- [Critical Analysis]
- [Extreme Analysis]
- [Trend Analysis]
- [Drought Analysis]
- [Conclusion]

1) Introduction

Climate change is a pressing global concern, impacting various environmental facets, including precipitation patterns. To comprehend the potential implications of climate change on precipitation, this comprehensive code conducts trend analyses, historical extreme precipitation assessments, and drought analyses. The code utilizes both historical and simulated future precipitation datasets to draw insights into changing precipitation trends, extreme events, and potential drought occurrences. Through visualizations and statistical models, the code aims to provide a robust understanding of how precipitation patterns have evolved historically and the potential trajectories they may follow in the future.
  
2) Reading netcdf function:

Reads a NetCDF file given the file path, latitude, and longitude of a city.
Extracts latitude, longitude, time, and precipitation data from the NetCDF file.
Identifies the closest grid point to the specified city.
Retrieves the precipitation time series for the closest point.
Converts precipitation units to daily totals (mm/day).
Returns a time series (xts object) with dates and corresponding daily precipitation values.

[process_directory function:]
Processes all NetCDF files in a given directory.
Applies process_netcdf to each file, resulting in a list of time series.
Merges the list of time series into a single time series using the merge function.
Returns the merged time series.

3) Dataset Analysis

[Data Preparation:]
The original precipitation data (gpr) is used, where each row represents a specific date along with its corresponding precipitation value.
The script converts the date information to the appropriate format using strptime and format functions.

[Aggregation:]
The aggregate function is employed to calculate the sum of daily precipitation for each unique date. This results in a new data frame (gpr.daily) where each row corresponds to a date, and the "pr" column contains the aggregated daily precipitation.

[Plotting:]
The script includes two plots to visualize the results:
The first plot displays the original daily precipitation data over time, providing insights into the day-to-day variability.
The second plot showcases the aggregated time series, presenting a smoothed representation of daily precipitation trends.

[Hydroplot Visualization:]
The hydroplot function is utilized to create a hydrograph-style plot for the aggregated time series. This plot offers a comprehensive view of the precipitation patterns, helping identify any notable trends or variations.

[Seasonal Analysis:]
The script performs seasonal analysis by categorizing precipitation into seasons (Winter, Spring, Summer, Autumn).
It creates a bar plot to visualize the seasonal variation in precipitation over the years.


4) Critical Analysis
This section of the code performs a critical analysis by comparing historical, observed, and future precipitation data. Here is a breakdown of the code:

[Monthly Summation:]
Calculates the total monthly precipitation for the historical and future datasets.
Reads observed precipitation data from a CSV file, processes it, and calculates the total monthly precipitation.

[Data Formatting:]
Formats historical, observed, and future monthly data into suitable structures.
Converts monthly data to a data frame with columns for Year, Month, and Precipitation.

[Yearly Aggregation:]
Aggregates monthly data into yearly totals for historical, observed, and future datasets.

[Plotting Yearly Totals:]
Creates a line plot comparing yearly precipitation totals for historical, observed, and future datasets.

[Monthly Data Comparison:]
Creates a line plot comparing historical, observed, and future monthly precipitation data over time.

[Scatter Plot and Evaluation Metrics:]
Generates a scatter plot comparing observed and predicted monthly precipitation.
Computes and displays Root Mean Squared Error (RMSE) and R-squared (R2) values for the scatter plot.

5) Extreme Analysis

[1-Day Maximum Rainfall:]

Calculates the maximum 1-day rainfall for each year in the future dataset.
Applies a linear regression model to explore trends in the annual maximum 1-day rainfall.
Creates a scatter plot with a fitted regression line to visualize the trend.

[7-Day Maximum Rainfall:]
Derives the 7-day period with the maximum precipitation for each year in the future dataset.
Performs linear regression on the total precipitation over the 7-day period to identify trends.
Presents a scatter plot with a fitted regression line for visualization.

[Historical Extreme Analysis (1-Day and 7-Day):]
Similar to the future analysis, explores historical data for both 1-day and 7-day extreme rainfall.
Provides scatter plots with linear regression lines for trend visualization.

[Exceedance Probability (Monthly Extreme):]
Ranks historical, observed, and future monthly precipitation values in descending order.
Computes exceedance probability based on the ranking.
Generates a plot with a log-scale x-axis to visualize the exceedance probability of extreme monthly precipitation.

6) Trend Analysis

[Annual Precipitation (Historical and Future):]
Computes the annual total precipitation for both historical and future datasets.
Creates scatter plots with dotted trend lines for both datasets.
Applies linear regression to quantify the trend in annual precipitation.
Annotates the plots with the slope of the regression line for historical and future periods.

[Future Time Windows (2005-2045, 2040-2070, 2060-2099):]
Divides the future dataset into three time windows.
Calculates the annual total precipitation for each time window.
Generates scatter plots with fitted regression lines and annotates the slope of each line.

[Historical Time Windows (1950-1980, 1970-2005):]
Similar to the future analysis, divides the historical dataset into two time windows.
Calculates the annual total precipitation for each time window.
Creates scatter plots with fitted regression lines and annotates the slope of each line.
   
7) Drought Analysis

[Future Drought:]
Extracts the year and month from the date in the future precipitation dataset.
Computes monthly total precipitation.
Calculates the Standardized Precipitation Index (SPI) using a 12-month aggregation period.
Creates a bar plot of SPI with different colors representing positive and negative values.
Includes horizontal dashed lines to indicate various drought severity levels (Moderate, Severe, Extreme).

[Historical Drought:]
Similar to the future analysis, processes the historical precipitation dataset.
Computes SPI and visualizes the historical precipitation deficit using a bar plot.
Includes the same drought severity level indicators.

[Common Features:]
Both analyses use the SPI, a widely used drought index.
The bar plots help visualize periods of drought and severity levels.
Annotations provide information about the severity of drought for specific SPI thresholds.

8) Conclusion

In conclusion, the code delivers a thorough examination of precipitation patterns, encompassing historical trends, extreme events, and potential future scenarios. The trend analyses offer insights into the historical changes in precipitation, revealing any discernible patterns or shifts. The assessment of extreme precipitation events, both historically and in the future, facilitates the understanding of the changing dynamics and potential risks associated with intense rainfall. Furthermore, the drought analyses, employing the Standardized Precipitation Index (SPI), contribute to assessing the vulnerability to drought conditions. By visualizing and categorizing drought severity levels, the code aids in identifying periods of concern and potential impacts on water resources and ecosystems.

This code serves as a valuable tool for climate scientists, environmental researchers, and policymakers seeking to understand and address the challenges posed by changing precipitation patterns. The visualizations and analyses provided empower stakeholders to make informed decisions regarding climate adaptation and mitigation strategies in the face of evolving precipitation dynamics.

