library(ncdf4)
library(sp)
library(geosphere)
library(xts)
library(tidyr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(hydroTSM)
library(forecast)
library(extRemes)

process_netcdf <- function(file_path, lat_city, lon_city) {
  ncin <- nc_open(file_path, write = FALSE, readunlim = TRUE, verbose = TRUE, auto_GMT = TRUE, suppress_dimvals = FALSE)
  
  lat <- ncvar_get(ncin, "lat")
  lon <- ncvar_get(ncin, "lon")
  
  mesh_size <- dim(lat)
  t <- ncvar_get(ncin, "time")
  obsdatadates <- as.Date(t, origin = '1949-12-01')
  
  coords_matrix <- cbind(as.vector(lon), as.vector(lat))
  mesh_points <- SpatialPoints(coords_matrix, proj4string = CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  city_point <- SpatialPoints(matrix(c(lon_city, lat_city), nrow = 1, ncol = 2), proj4string = CRS("+proj=longlat +datum=WGS84"), bbox = NULL)
  
  closest_point <- which.min(distGeo(mesh_points, city_point))
  point_ID_matrix <- matrix(1:length(lat), nrow = nrow(lat), ncol = ncol(lat))
  ij <- which(point_ID_matrix == closest_point, arr.ind = TRUE)
  timeseries <- ncvar_get(ncin, varid = 'pr', start = c(ij[1], ij[2], 1), count = c(1, 1, -1))
  
  
  timeseries <- timeseries * (86400 * 1)
  
  nc_close(ncin)
  
  return(xts(timeseries, obsdatadates))
}

process_directory <- function(directory_path, lat_city, lon_city) {
  setwd(directory_path)
  nc_files <- list.files(pattern = "\\.nc$")
  timeseries_list <- lapply(nc_files, function(file) process_netcdf(file, lat_city, lon_city))
  return(do.call(merge, timeseries_list))
}

process_precipitation <- function(precipitation_df) {
  precipitation_df <- data.frame(Date = index(precipitation_df), coredata(precipitation_df))
  precipitation_df[] <- lapply(precipitation_df, as.character)
  precipitation_df$pr <- apply(precipitation_df[, -1, drop = FALSE], 1, function(x) paste(na.omit(x), collapse = ""))
  return(precipitation_df[, c("Date", "pr"), drop = FALSE])
}

calculate_monthly_avg <- function(precipitation_df) {
  precipitation_df <- precipitation_df %>%
    mutate(Month = format(Date, "%Y-%m"))
  
  monthly_sum <- precipitation_df %>%
    group_by(Month) %>%
    summarise(Total_Precipitation = sum(pr))
  
  monthly_sum <- monthly_sum %>%
    mutate(Month = substr(Month, 6, 7))
  
  monthly_avg <- monthly_sum %>%
    group_by(Month) %>%
    summarise(Average_Precipitation = sum(Total_Precipitation) / length(unique(lubridate::year(precipitation_df$Date))))
  
  
  monthly_avg$Month <- month.name[as.numeric(monthly_avg$Month)]
  
  return(monthly_avg)
}


Historical_Precipitation_df <- process_directory("/Users/cruel-mac/Downloads/Drought-assignment/Southasia/historical_chennai/precipitation/", 13.0827, 80.2707)
Historical_Precipitation_df <- process_precipitation(Historical_Precipitation_df)

Future_Precipitation_df <- process_directory("/Users/cruel-mac/Downloads/Drought-assignment/Southasia/future_RCP-8.5_south/Precip_RCP8.5/", 13.0827, 80.2707)
Future_Precipitation_df <- process_precipitation(Future_Precipitation_df)


Historical_Precipitation_df$Date <- as.Date(Historical_Precipitation_df$Date)
Historical_Precipitation_df$pr <- as.numeric(Historical_Precipitation_df$pr)
Future_Precipitation_df$Date <- as.Date(Future_Precipitation_df$Date)
Future_Precipitation_df$pr <- as.numeric(Future_Precipitation_df$pr)



Historical_monthly_avg <- calculate_monthly_avg(Historical_Precipitation_df)
Future_monthly_avg<-calculate_monthly_avg(Future_Precipitation_df)


#---------------------# TOTAL HISTORICAL DATASET ANALYSIS#---------------------#
gpr<-Historical_Precipitation_df


Date.gpr<-strptime(gpr$Date, format = "%Y-%m-%d ")

Dates.gpr<-format(Date.gpr, "%Y-%m-%d")

gpr.daily<-aggregate(gpr$pr, by=list(Dates.gpr), FUN=sum)
names(gpr.daily)

names(gpr.daily)<-c("Dates.gpr","pr")
gpr.daily$Dates.gpr=as.Date(gpr.daily$Dates.gpr,"%Y-%m-%d")
plot(Date.gpr,gpr$pr,type = "l", xlab="Year",ylab="Precipitation")
gpr.daily.ts=zoo(gpr.daily$pr,order.by = gpr.daily$Dates.gpr)
gpr.daily.ts

plot(gpr.daily.ts, xlab = "Year", ylab =" Precipitation")


hydroplot(gpr.daily.ts, var.type = "Precipitation", var.unit = "(mm)",xlab = "Time", ylab = "precipitation" )

# seasonal analysis


gpr$Date <- as.Date(gpr$Date)
gpr <- gpr %>%
  mutate(Year = lubridate::year(Date),
         Month = lubridate::month(Date, label = TRUE))
monthly_summary <- gpr %>%
  group_by(Year, Month) %>%
  summarise(Total_Precipitation = sum(pr))
correct_month_order <- month.abb[c(12, 1:11)]
monthly_summary$Month <- factor(monthly_summary$Month, levels = correct_month_order, ordered = TRUE)
monthly_summary_wide <- monthly_summary %>%
  pivot_wider(names_from = Month, values_from = Total_Precipitation, names_prefix = "Month_")


colnames(monthly_summary_wide) <- c("Year", correct_month_order)
monthly_summary_wide <- monthly_summary_wide[, c("Year", month.abb[c(1:11, 12)])]
monthly_summary_wide <- monthly_summary_wide %>%
  mutate(
    Winter = Jan + Feb + Dec,           # Winter: Dec + Jan + Feb
    Spring = Mar + Apr + May,           # Spring: Mar + Apr + May
    Summer = Jun + Jul + Aug,           # Summer: Jun + Jul + Aug
    Autumn = Sep + Oct + Nov            # Autumn: Sep + Oct + Nov
  )

plot_data <- tidyr::gather(monthly_summary_wide, Season, Value, Winter:Autumn)


ggplot(plot_data, aes(x = Year, y = Value, fill = Season)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Seasonal analysis - Historical",
       x = "Year",
       y = "Precipitation (mm)",
       fill = "Season") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12), 
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )


# comparing -historical vs future vs observed data


Historical_monthly_avg$Month <- factor(Historical_monthly_avg$Month, levels = month.name)
Future_monthly_avg$Month <- factor(Future_monthly_avg$Month, levels = month.name)


Real_observed_data <- read.csv("/Users/cruel-mac/Downloads/Drought-assignment/Real_data/csvFile94.csv", header = TRUE)

Real_observed_data$Dates <- as.Date(Real_observed_data$Dates)

Real_observed_data <- Real_observed_data %>%
  mutate(
    Year = year(Dates),
    Month = month(Dates)
  )


month_names <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")


Real_observed_data$Month <- factor(month(Real_observed_data$Month), levels = 1:12, labels = month_names)

Real_monthly_sum <- Real_observed_data %>%
  group_by(Month, Year) %>%
  summarise(Total_Precipitation = sum(ACTUAL..mm., na.rm = TRUE))

Real_monthly_avg <- Real_observed_data %>%
  group_by(Month) %>%
  summarise(Average_Precipitation = mean(ACTUAL..mm., na.rm = TRUE))
str(Real_monthly_avg)


combined_data <- bind_rows(
  mutate(Historical_monthly_avg, dataset = "Historical-(1950-2005)"),
  mutate(Future_monthly_avg, dataset = "Future-RCP8.5-( 2006-2098)"),
  mutate(Real_monthly_avg, dataset = "Observed-(1950-2023)")
)

ggplot(combined_data, aes(x = Month, y = Average_Precipitation, fill = dataset)) +
  geom_ribbon(aes(ymin = 0, ymax = Average_Precipitation), alpha = 0.5) +
  geom_line(aes(group = dataset), size = 1, linetype = "dashed") +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
  labs(title = "Monthly Mean Rainfall Distribution",
       x = "Month",
       y = "Precipitation (mm)",
       fill = "Dataset") +
  theme_minimal() +
  facet_wrap(~dataset, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("red", "blue", "black")) +  
  theme(
    axis.text = element_text(size = 10, color = "black"),  
    axis.title = element_text(size = 14, color = "black"),  
    axis.line = element_line(color = "black"),  
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
    legend.text = element_text(size = 11),  
    legend.title = element_text(size = 12)  
  )




#-----------------# CRITICAL ANALYSIS#------------------------------------------#

Historical_monthly_sum <- Historical_Precipitation_df %>%
  mutate(Month = format(Date, "%Y-%m")) %>%
  group_by(Month) %>%
  summarise(Total_Precipitation = sum(pr))


Future_monthly_sum <- Future_Precipitation_df %>%
  mutate(Month = format(Date, "%Y-%m")) %>%
  group_by(Month) %>%
  summarise(Total_Precipitation = sum(pr))


Real_observed_data <- read.csv("/Users/cruel-mac/Downloads/Drought-assignment/Real_data/csvFile94.csv", header = TRUE)

Real_observed_data$Dates <- as.Date(Real_observed_data$Dates)

Real_observed_data <- Real_observed_data %>%
  mutate(
    Year = year(Dates),
    Month = month(Dates)
  )


Real_monthly_data <- Real_observed_data[, c("Year", "Month", "ACTUAL..mm.")]

Real_monthly_data <- Real_monthly_data %>%
  rename(Actual_pr = ACTUAL..mm.)

historical_monthly_data <- data.frame(
  Year = as.numeric(substr(Historical_monthly_sum$Month, 1, 4)),
  Month = as.numeric(substr(Historical_monthly_sum$Month, 6, 7)),
  Precipitation = Historical_monthly_sum$Total_Precipitation
)


future_monthly_data <- data.frame(
  Year = as.numeric(substr(Future_monthly_sum$Month, 1, 4)),
  Month = as.numeric(substr(Future_monthly_sum$Month, 6, 7)),
  Precipitation = Future_monthly_sum$Total_Precipitation
)


str(historical_monthly_data)
str(future_monthly_data)
str(Real_monthly_data)



historical_monthly_data$Date <- as.Date(paste(historical_monthly_data$Year, historical_monthly_data$Month, 1, sep = "-"), format = "%Y-%m-%d")
future_monthly_data$Date <- as.Date(paste(future_monthly_data$Year, future_monthly_data$Month, 1, sep = "-"), format = "%Y-%m-%d")
Real_monthly_data$Date <- as.Date(paste(Real_monthly_data$Year, Real_monthly_data$Month, 1, sep = "-"), format = "%Y-%m-%d")

historical_yearly_data <- historical_monthly_data %>%
  group_by(Year) %>%
  summarise(Total_Precipitation = sum(Precipitation))

future_yearly_data <- future_monthly_data %>%
  group_by(Year) %>%
  summarise(Total_Precipitation = sum(Precipitation))

Real_yearly_data <- Real_monthly_data %>%
  group_by(Year) %>%
  summarise(Total_Precipitation = sum(Actual_pr))

str(historical_yearly_data)
str(future_yearly_data)
str(Real_yearly_data)

future_yearly_data_subset <- subset(future_yearly_data, Year <= 2023)


ggplot() +
  geom_line(data = historical_yearly_data, aes(x = Year, y = Total_Precipitation, color = "Historical"), size = 1.5) +
  geom_line(data = Real_yearly_data, aes(x = Year, y = Total_Precipitation, color = "Real"), size = 1.5) +
  geom_line(data = future_yearly_data_subset, aes(x = Year, y = Total_Precipitation, color = "Future"), size = 1.5) +
  labs(title = "Historical vs. Real vs. Future Yearly Precipitation",
       x = "Year",
       y = "Total Precipitation") +
  scale_color_manual(values = c("Historical" = "blue", "Real" = "black", "Future" = "red")) +
  theme_minimal()

future_monthly_data_subset <- subset(future_monthly_data, Year <= 2023)

ggplot() +
  geom_line(data = historical_monthly_data, aes(x = Date, y = Precipitation, color = "Historical"), size = 0.3, alpha = 0.8) +
  geom_line(data = Real_monthly_data, aes(x = as.Date(paste(Year, Month, 1, sep = "-"), format = "%Y-%m-%d"), y = Actual_pr, color = "Observed"), size = 0.3, alpha = 0.8) +
  geom_line(data = future_monthly_data_subset, aes(x = as.Date(paste(Year, Month, 1, sep = "-"), format = "%Y-%m-%d"), y = Precipitation, color = "RCP-8.5"), size = 0.3, alpha = 0.8) +
  labs(title = "Observed vs Climate Model", face = "bold",
       x = "Date",
       y = "Precipitation (mm/month)") +
  scale_color_manual(values = c("Historical" = "blue", "Observed" = "black", "RCP-8.5" = "red"),
                     breaks = c("Historical", "Observed", "RCP-8.5")) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12),  
    plot.title = element_text(size = 16)  
  )


scatter_data <- merge(historical_monthly_data, Real_monthly_data, by = c("Year", "Month"))



scatter_data_filtered <- subset(scatter_data, Actual_pr <= 1500 & Precipitation <= 1500)


rmse_value <- sqrt(mean((scatter_data_filtered$Actual_pr - scatter_data_filtered$Precipitation)^2))
r2_value <- cor(scatter_data_filtered$Actual_pr, scatter_data_filtered$Precipitation)^2


ggplot(scatter_data_filtered, aes(x = Actual_pr, y = Precipitation)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "gray") +  
  labs(
    title = sprintf("1950-2023\nRMSE: %.2f, R2: %.2f", rmse_value, r2_value),
    x = "Observed Precipitation (mm/month)",
    y = "Predicted Precipitation (mm/month)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),  
    plot.title = element_text(size = 16),  
    axis.line = element_line(color = "black", size = 1.5)  
  ) +
  coord_fixed(ratio = 1) +  
  xlim(0, 1500) +  
  ylim(0, 1500)  


str(historical_monthly_data)
str(Real_monthly_data)
str(future_monthly_data)


#-------------------------# EXTREME ANALYSIS#------------------------------------#


#future_extreme
Future_Precipitation_df$Year <- lubridate::year(Future_Precipitation_df$Date)

max_rainfall_per_year <- Future_Precipitation_df %>%
  group_by(Year) %>%
  summarise(Max_Rainfall = max(pr),
            Date_of_Max_Rainfall = Date[which.max(pr)])
linear_model <- lm(Max_Rainfall ~ Date_of_Max_Rainfall, data = max_rainfall_per_year)
slope_value <- coef(linear_model)[2]
scatter_plot <- ggplot(max_rainfall_per_year, aes(x = Date_of_Max_Rainfall, y = Max_Rainfall)) +
  geom_point(color = "blue", size = 2, alpha = 1.3) +  
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dotdash") +  
  geom_text(x = Inf, y = Inf, label = paste("Slope: ", round(slope_value, 4)), 
            hjust = 1, vjust = 1, size = 4, color = "black") +  
  labs(title = "1-Day Maximum Rainfall",
       x = "Date",
       y = "Precipitation (mm)",
       color = "Legend") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 11, color = "black"),  
    axis.line = element_line(color = "black")  
  )


print(scatter_plot)


max_1day_rainfall_per_year <- Future_Precipitation_df %>%
  group_by(Year) %>%
  summarise(Max_1Day_Rainfall = max(pr),
            Date_of_Max_1Day_Rainfall = Date[which.max(pr)])


max_7day_rainfall_per_year <- max_1day_rainfall_per_year %>%
  mutate(
    Date_start = Date_of_Max_1Day_Rainfall - days(3),
    Date_end = Date_of_Max_1Day_Rainfall + days(3)
  ) %>%
  left_join(Future_Precipitation_df, by = "Year") %>%
  filter(Date >= Date_start & Date <= Date_end) %>%
  group_by(Year) %>%
  summarise(Total_7Day_Rainfall = sum(pr))


print(max_7day_rainfall_per_year)


linear_model_7day <- lm(Total_7Day_Rainfall ~ Year, data = max_7day_rainfall_per_year)


slope_value_7day <- coef(linear_model_7day)[2]


scatter_plot_7day <- ggplot(max_7day_rainfall_per_year, aes(x = Year, y = Total_7Day_Rainfall)) +
  geom_point(color = "purple", size = 2, alpha = 1.3) +  
  geom_smooth(method = "lm", se = FALSE, color = "brown", linetype = "dotdash") +  
  geom_text(x = Inf, y = Inf, label = paste("Slope: ", round(slope_value_7day, 4)), 
            hjust = 1, vjust = 1, size = 4, color = "black") + 
  labs(title = "7-Day maximum Rainfall",
       x = "Year",
       y = "Precipitation (mm)",
       color = "Legend") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 11, color = "black"),  
    axis.line = element_line(color = "black")  
  )

print(scatter_plot_7day)


#historical extreme

Historical_Precipitation_df$Year <- lubridate::year(Historical_Precipitation_df$Date)

Historical_max_rainfall_per_year <- Historical_Precipitation_df %>%
  group_by(Year) %>%
  summarise(Max_Rainfall = max(pr),
            Date_of_Max_Rainfall = Date[which.max(pr)])
print(Historical_max_rainfall_per_year)
linear_model <- lm(Max_Rainfall ~ Date_of_Max_Rainfall, data = Historical_max_rainfall_per_year)
slope_value <- coef(linear_model)[2]


scatter_plot <- ggplot(Historical_max_rainfall_per_year, aes(x = Date_of_Max_Rainfall, y = Max_Rainfall)) +
  geom_point(color = "blue", size = 2, alpha = 1.3) +  
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dotdash") +  
  geom_text(x = Inf, y = Inf, label = paste("Slope: ", round(slope_value, 4)), 
            hjust = 1, vjust = 1, size = 4, color = "black") +  
  labs(title = "1-Day Maximum Rainfall",
       x = "Date",
       y = "Precipitation (mm)",
       color = "Legend") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 11, color = "black"),  
    axis.line = element_line(color = "black")  
  )


print(scatter_plot)

max_1day_rainfall_per_year_historical <- Historical_Precipitation_df %>%
  group_by(Year) %>%
  summarise(Max_1Day_Rainfall = max(pr),
            Date_of_Max_1Day_Rainfall = Date[which.max(pr)])


max_7day_rainfall_per_year_historical <- max_1day_rainfall_per_year_historical %>%
  mutate(
    Date_start = Date_of_Max_1Day_Rainfall - days(3),
    Date_end = Date_of_Max_1Day_Rainfall + days(3)
  ) %>%
  left_join(Historical_Precipitation_df, by = "Year") %>%
  filter(Date >= Date_start & Date <= Date_end) %>%
  group_by(Year) %>%
  summarise(Total_7Day_Rainfall = sum(pr))


linear_model_7day_historical <- lm(Total_7Day_Rainfall ~ Year, data = max_7day_rainfall_per_year_historical)


slope_value_7day_historical <- coef(linear_model_7day_historical)[2]


scatter_plot_7day_historical <- ggplot(max_7day_rainfall_per_year_historical, aes(x = Year, y = Total_7Day_Rainfall)) +
  geom_point(color = "purple", size = 2, alpha = 1.3) +  
  geom_smooth(method = "lm", se = FALSE, color = "brown", linetype = "dotdash") +  
  geom_text(x = Inf, y = Inf, label = paste("Slope: ", round(slope_value_7day_historical, 4)), 
            hjust = 1, vjust = 1, size = 4, color = "black") +  
  labs(title = "7-Day Maximum Rainfall",
       x = "Year",
       y = "Precipitation (mm)",
       color = "Legend") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 13, hjust = 0.5),  
    axis.title = element_text(size = 12),  
    axis.text = element_text(size = 11, color = "black"),  
    axis.line = element_line(color = "black")  
  )


print(scatter_plot_7day_historical)

# Exceedance probability- For Monthly extreme
historical_monthly_data_prob <- historical_monthly_data %>%
  arrange(desc(Precipitation)) %>%
  mutate(Exceedance_Probability = (1:length(Precipitation)) / length(Precipitation))

historical_monthly_data_prob <- historical_monthly_data %>%
  arrange(desc(Precipitation)) %>%
  mutate(Exceedance_Probability = (1:length(Precipitation)) / length(Precipitation))

real_monthly_data_prob <- Real_monthly_data %>%
  arrange(desc(Actual_pr)) %>%
  mutate(Exceedance_Probability = (1:length(Actual_pr)) / length(Actual_pr))

future_monthly_data_prob <- future_monthly_data %>%
  arrange(desc(Precipitation)) %>%
  mutate(Exceedance_Probability = (1:length(Precipitation)) / length(Precipitation))


historical_monthly_data_prob <- rename(historical_monthly_data_prob, Precipitation = Precipitation)
real_monthly_data_prob <- rename(real_monthly_data_prob, Precipitation = Actual_pr)
future_monthly_data_prob <- rename(future_monthly_data_prob, Precipitation = Precipitation)


historical_monthly_data_prob <- select(historical_monthly_data_prob, Exceedance_Probability, Precipitation)
real_monthly_data_prob <- select(real_monthly_data_prob, Exceedance_Probability, Precipitation)
future_monthly_data_prob <- select(future_monthly_data_prob, Exceedance_Probability, Precipitation)



combined_data <- bind_rows(
  mutate(historical_monthly_data_prob, Source = "Historical"),
  mutate(real_monthly_data_prob, Source = "Real"),
  mutate(future_monthly_data_prob, Source = "Future")
)


str(combined_data)

summary(real_monthly_data_prob)
summary(future_monthly_data_prob)


ggplot(combined_data, aes(x = Exceedance_Probability, y = Precipitation, color = Source)) +
  geom_point() +
  geom_line(aes(group = Source), linetype = "solid", size = 1) +  
  labs(title = "Exceedance Probability- Monthly Extreme",
       x = "Exceedance Probability (Log10)",
       y = "Precipitation (mm)") +
  theme_minimal() +
  scale_x_log10() +  # Use log scale for x-axis
  scale_y_continuous(breaks = seq(0, max(combined_data$Precipitation), by = 200)) +  
  scale_color_manual(values = c("red", "blue", "black"), labels = c("Future", "Historical", "Observed")) + 
  theme(
    axis.line = element_line(color = "black", size = 1.5),  
    axis.text = element_text(size = 12),  
    plot.title = element_text(size = 14, hjust = 0.5)  
  )


#--------------------------------# TREND ANALYSIS#------------------------------#

historical_annual_total <- Historical_Precipitation_df%>%
  mutate(Year = lubridate::year(Date)) %>%
  group_by(Year) %>%
  summarise(Total_Precipitation = sum(pr, na.rm = TRUE))


future_annual_total <- Future_Precipitation_df %>%
  mutate(Year = lubridate::year(Date)) %>%
  group_by(Year) %>%
  summarise(Total_Precipitation = sum(pr, na.rm = TRUE))

ggplot() +
  geom_point(data = historical_annual_total, aes(x = Year, y = Total_Precipitation), color = "blue") +
  geom_line(data = historical_annual_total, aes(x = Year, y = Total_Precipitation), linetype = "dotted", color = "blue") +
  geom_point(data = future_annual_total, aes(x = Year, y = Total_Precipitation), color = "red") +
  geom_line(data = future_annual_total, aes(x = Year, y = Total_Precipitation), linetype = "dotted", color = "red") +
  geom_smooth(data = historical_annual_total, aes(x = Year, y = Total_Precipitation), method = "lm", se = FALSE, color = "blue") +
  geom_smooth(data = future_annual_total, aes(x = Year, y = Total_Precipitation), method = "lm", se = FALSE, color = "red") +
  labs(title = "Annual Precipitation (Historical and Future)",
       x = "Year",
       y = "Precipitation (mm)") +
  theme_minimal() +
  annotate("text", x = 1950, y = max(historical_annual_total$Total_Precipitation), 
           label = paste("Slope (Historical):", round(coef(lm(Total_Precipitation ~ Year, data = historical_annual_total))[2], 4), "mm per year"), 
           color = "blue", hjust = 0) +
  annotate("text", x = 2021, y = max(future_annual_total$Total_Precipitation), 
           label = paste("Slope (Future):", round(coef(lm(Total_Precipitation ~ Year, data = future_annual_total))[2], 4), "mm per year"), 
           color = "red", hjust = 0) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),)

#Future Time windows

future_precipitation_window1 <-Future_Precipitation_df %>%
  filter(Date >= "2005-01-01" & Date <= "2045-12-31") %>%
  group_by(Year = lubridate::year(Date)) %>%
  summarise(Total_Precipitation = sum(pr, na.rm = TRUE))

future_precipitation_window2 <- Future_Precipitation_df %>%
  filter(Date >= "2040-01-01" & Date <= "2070-12-31") %>%
  group_by(Year = lubridate::year(Date)) %>%
  summarise(Total_Precipitation = sum(pr, na.rm = TRUE))

future_precipitation_window3 <- Future_Precipitation_df %>%
  filter(Date >= "2060-01-01" & Date <= "2099-12-31") %>%
  group_by(Year = lubridate::year(Date)) %>%
  summarise(Total_Precipitation = sum(pr, na.rm = TRUE))

ggplot() +
  geom_point(data = future_precipitation_window1, aes(x = Year, y = Total_Precipitation), color = "red", alpha = 0.7) +
  geom_line(data = future_precipitation_window1, aes(x = Year, y = Total_Precipitation), linetype = "solid", color = "red", alpha = 0.7) +
  geom_smooth(data = future_precipitation_window1, aes(x = Year, y = Total_Precipitation), method = "lm", se = FALSE, linetype = "dotted", color = "red") +
  labs(title = "Year (2005-2045)",
       x = "Year",
       y = " Annual Precipitation (mm)") +
  theme_minimal() +
  annotate("text", x = 2005, y = max(future_precipitation_window1$Total_Precipitation), 
           label = paste(" y =", round(coef(lm(Total_Precipitation ~ Year, data = future_precipitation_window1))[2], 4), "x +", round(coef(lm(Total_Precipitation ~ Year, data = future_precipitation_window1))[1], 4)), 
           color = "black", size = 6, hjust = 0, vjust = 0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))

ggplot() +
  geom_point(data = future_precipitation_window2, aes(x = Year, y = Total_Precipitation), color = "red", alpha = 0.7) +
  geom_line(data = future_precipitation_window2, aes(x = Year, y = Total_Precipitation), linetype = "solid", color = "red", alpha = 0.7) +
  geom_smooth(data = future_precipitation_window2, aes(x = Year, y = Total_Precipitation), method = "lm", se = FALSE, linetype = "dotted", color = "red") +
  labs(title = "Year (2040-2070)",
       x = "Year",
       y = "Annual Precipitation (mm)") +
  theme_minimal() +
  annotate("text", x = 2042, y = max(future_precipitation_window2$Total_Precipitation), 
           label = paste(" y =", round(coef(lm(Total_Precipitation ~ Year, data = future_precipitation_window2))[2], 4), "x +", round(coef(lm(Total_Precipitation ~ Year, data = future_precipitation_window2))[1], 4)), 
           color = "black", size = 6, hjust = 0, vjust = 0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))

ggplot() +
  geom_point(data = future_precipitation_window3, aes(x = Year, y = Total_Precipitation), color = "red", alpha = 0.7) +
  geom_line(data = future_precipitation_window3, aes(x = Year, y = Total_Precipitation), linetype = "solid", color = "red", alpha = 0.7) +
  geom_smooth(data = future_precipitation_window3, aes(x = Year, y = Total_Precipitation), method = "lm", se = FALSE, linetype = "dotted", color = "red") +
  labs(title = " Year(2060-2099)",
       x = "Year",
       y = "Annual Precipitation (mm)") +
  theme_minimal() +
  annotate("text", x = 2063, y = max(future_precipitation_window3$Total_Precipitation), 
           label = paste(" y =", round(coef(lm(Total_Precipitation ~ Year, data = future_precipitation_window3))[2], 4), "x +", round(coef(lm(Total_Precipitation ~ Year, data = future_precipitation_window3))[1], 4)), 
           color = "black", size = 6, hjust = 0, vjust = 0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))



# historical time windows

historical_precipitation_window1 <- Historical_Precipitation_df %>%
  filter(Date >= "1950-01-01" & Date <= "1980-12-31") %>%
  group_by(Year = lubridate::year(Date)) %>%
  summarise(Total_Precipitation = sum(pr, na.rm = TRUE))

historical_precipitation_window2 <- Historical_Precipitation_df %>%
  filter(Date >= "1970-01-01" & Date <= "2005-12-31") %>%
  group_by(Year = lubridate::year(Date)) %>%
  summarise(Total_Precipitation = sum(pr, na.rm = TRUE))


ggplot() +
  geom_point(data = historical_precipitation_window1, aes(x = Year, y = Total_Precipitation), color = "blue", alpha = 0.7) +
  geom_line(data = historical_precipitation_window1, aes(x = Year, y = Total_Precipitation), linetype = "solid", color = "blue", alpha = 0.7) +
  geom_smooth(data = historical_precipitation_window1, aes(x = Year, y = Total_Precipitation), method = "lm", se = FALSE, linetype = "dotted", color = "blue") +
  labs(title = "Year (1950-1980)",
       x = "Year",
       y = " Annual Precipitation (mm)") +
  theme_minimal() +
  annotate("text", x = 1950, y = max(historical_precipitation_window1$Total_Precipitation), 
           label = paste(" y =", round(coef(lm(Total_Precipitation ~ Year, data = historical_precipitation_window1))[2], 4), "x +", round(coef(lm(Total_Precipitation ~ Year, data = historical_precipitation_window1))[1], 4)), 
           color = "black", size = 6, hjust = 0, vjust = 0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))


ggplot() +
  geom_point(data = historical_precipitation_window2, aes(x = Year, y = Total_Precipitation), color = "blue", alpha = 0.7) +
  geom_line(data = historical_precipitation_window2, aes(x = Year, y = Total_Precipitation), linetype = "solid", color = "blue", alpha = 0.7) +
  geom_smooth(data = historical_precipitation_window2, aes(x = Year, y = Total_Precipitation), method = "lm", se = FALSE, linetype = "dotted", color = "blue") +
  labs(title = " Year (1970-2005)",
       x = "Year",
       y = " Annual Precipitation (mm)") +
  theme_minimal() +
  annotate("text", x = 1970, y = max(historical_precipitation_window2$Total_Precipitation), 
           label = paste(" y =", round(coef(lm(Total_Precipitation ~ Year, data = historical_precipitation_window2))[2], 4), "x +", round(coef(lm(Total_Precipitation ~ Year, data = historical_precipitation_window2))[1], 4)), 
           color = "black", size = 6, hjust = 0, vjust = 0.2) +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))





#-----------------------------#DROUGHT ANALYSIS#-------------------------------#

#Future_drought

Future_Precipitation_df <- Future_Precipitation_df %>%
  mutate(Year = format(Date, "%Y"),
         Month = format(Date, "%m"))

monthly_precipitation <- Future_Precipitation_df %>%
  group_by(Year, Month) %>%
  summarise(Total_Precipitation = sum(pr))


print(monthly_precipitation)
SPI<-spi(monthly_precipitation$Total_Precipitation,12)
df <- fortify.zoo(SPI$fitted)
names(df) <- c("ID", "SPI")
df <- df %>%
  select(-ID) %>%
  mutate(Period = as.yearmon(paste(monthly_precipitation$Year, monthly_precipitation$Month, "01", sep = "-")),
         sign = ifelse(SPI > 0, "pos", "neg")) %>%
  na.omit()

head(df, 12)

str(df)

ggplot(df) +
  geom_bar(aes(x = Period, y = SPI, fill = sign), stat = "identity") +
  scale_fill_manual(values = c("pos" = "darkblue", "neg" = "red")) +
  scale_y_continuous(limits = c(-3.7, 3.5), breaks = -3:3.5) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(2005, 2099, 5)) +
  labs(y = "SPI 12", x = "Year") +
  ggtitle("Future Precipitation Deficit") +
  theme(axis.title.x = element_text(face = "bold.italic", size = 11),
        axis.title.y = element_text(face = "bold.italic", size = 11),
        axis.text.x = element_text(size = 11, color = "#000000", hjust = 0.6),
        axis.text.y = element_text(size = 11, color = "#000000"),
        axis.line = element_line(size = 0.5, color = "#000000"),
        plot.title = element_text(vjust = 1.5, hjust = 0.5, size = 14, face = "bold"),  
        legend.position = "none") +
  geom_hline(yintercept = c(-1, -1.5, -2, -2.5),
             linetype = "dashed",
             color = "gray") +
  geom_text(aes(x = max(df$Period) + 0.1, y = -1.25, label = "Moderate Drought"), hjust = 1.1, vjust = 0, color = "black") +
  geom_text(aes(x = max(df$Period) + 0.1, y = -1.75, label = "Severe Drought"), hjust = 1.1, vjust = 0, color = "black") +
  geom_text(aes(x = max(df$Period) + 0.1, y = -2.25, label = "Extreme Drought"), hjust = 1.1, vjust = 0, color = "black")



#Historical Drought

Historical_Precipitation_df_1 <- Historical_Precipitation_df %>%
  mutate(Year = format(Date, "%Y"),
         Month = format(Date, "%m"))


monthly_precipitation <- Historical_Precipitation_df_1 %>%
  group_by(Year, Month) %>%
  summarise(Total_Precipitation = sum(pr))

print(monthly_precipitation)

SPI<-spi(monthly_precipitation$Total_Precipitation,12)


df <- fortify.zoo(SPI$fitted)


names(df) <- c("ID", "SPI")


df <- df %>%
  select(-ID) %>%
  mutate(Period = as.yearmon(paste(monthly_precipitation$Year, monthly_precipitation$Month, "01", sep = "-")),
         sign = ifelse(SPI > 0, "pos", "neg")) %>%
  na.omit()

head(df, 12)

str(df)


ggplot(df) +
  geom_bar(aes(x = Period, y = SPI, fill = sign), stat = "identity") +
  scale_fill_manual(values = c("pos" = "darkblue", "neg" = "red")) +
  scale_y_continuous(limits = c(-3.7, 3.5), breaks = -3:3.5) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(1950, 2005, 5)) +
  labs(y = "SPI 12", x = "Year") +
  ggtitle("Historical Precipitation Deficit") +
  theme(axis.title.x = element_text(face = "bold.italic", size = 11),
        axis.title.y = element_text(face = "bold.italic", size = 11),
        axis.text.x = element_text(size = 11, color = "#000000", hjust = 0.6),
        axis.text.y = element_text(size = 11, color = "#000000"),
        axis.line = element_line(size = 0.5, color = "#000000"),
        plot.title = element_text(vjust = 1.5, hjust = 0.5, size = 14, face = "bold"),  
        legend.position = "none") +
  geom_hline(yintercept = c(-1, -1.5, -2, -2.5),
             linetype = "dashed",
             color = "gray") +
  geom_text(aes(x = max(df$Period) + 0.1, y = -1.25, label = "Moderate Drought"), hjust = 1.1, vjust = 0, color = "black") +
  geom_text(aes(x = max(df$Period) + 0.1, y = -1.75, label = "Severe Drought"), hjust = 1.1, vjust = 0, color = "black") +
  geom_text(aes(x = max(df$Period) + 0.1, y = -2.25, label = "Extreme Drought"), hjust = 1.1, vjust = 0, color = "black")




summary(SPI)
names(SPI)
SPI$call
SPI$fitted
SPI$coefficients
