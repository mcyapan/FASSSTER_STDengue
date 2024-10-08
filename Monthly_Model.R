library(dplyr)
library(INLA)
library(sf)
library(ggplot2)
library(tidyr)
library(dlnm)

############################################################################### load data 

load("load.RData")
load("districts.rdata")
load("monthly_dengue_data.rdata")

monthly_dengue_data <- read.csv("monthly_dengue_data.csv")

View(monthly_dengue_data)
# data is complete until 7/2024
# model fitting will be done until 6/2024, then prediction will be done on 7/2024

############################################################################### cleaning data format

# creating lags in dataset

monthly_dengue_data <- monthly_dengue_data %>%
  group_by(District) %>%
  mutate(lag_precip_2 = lag(total_precip, 2),
         lag_tmax_2 = lag(avg_tempmax, 2), 
         lag_tmin_2=lag(avg_tempmin,2), 
         lag_hum_2=lag(avg_humidity,2)) %>%
  ungroup()

monthly_dengue_data <- monthly_dengue_data %>%
  group_by(District) %>%
  mutate(total_cases_lag_1 = lag(total_cases, 1),
         popden_lag_1 = lag(PopDen, 1),
         hh_lag_1 = lag(Average_Household, 1)) %>%
  ungroup()

# remove data from 8/2024 onwards

monthly_dengue_data <- monthly_dengue_data %>%
  filter(year < 2024 | (year == 2024 & month <8))

# setting 7/2024 data to NA (for prediction)

monthly_dengue_data_july = monthly_dengue_data$total_cases[monthly_dengue_data$year == 2024 & monthly_dengue_data$month == 7]
monthly_dengue_data_nojuly = monthly_dengue_data
monthly_dengue_data_nojuly$total_cases[monthly_dengue_data_nojuly$year == 2024 & monthly_dengue_data_nojuly$month == 7] <- NA

View(monthly_dengue_data_july)
View(monthly_dengue_data_nojuly)


############################################################################### constructing cross basis matrix

# set minimum lag to 3 months and maximum lag to 12 months

minlag = 3
maxlag = 12

# construct lag vectors

lag_precip = tsModel::Lag(monthly_dengue_data_nojuly$total_precip, 
                          group = monthly_dengue_data_nojuly$District, 
                          k = minlag:maxlag)
lag_tmax = tsModel::Lag(monthly_dengue_data_nojuly$avg_tempmax, 
                        group = monthly_dengue_data_nojuly$District, 
                        k = minlag:maxlag)
lag_tmin = tsModel::Lag(monthly_dengue_data_nojuly$avg_tempmin, 
                        group = monthly_dengue_data_nojuly$District, 
                        k = minlag:maxlag)
lag_hum = tsModel::Lag(monthly_dengue_data_nojuly$avg_humidity, 
                       group = monthly_dengue_data_nojuly$District, 
                       k = minlag:maxlag)
lag_tcases = tsModel::Lag(monthly_dengue_data_nojuly$total_cases, 
                          group = monthly_dengue_data_nojuly$District, 
                          k = minlag:maxlag)
lag_popden = tsModel::Lag(monthly_dengue_data_nojuly$PopDen, 
                          group = monthly_dengue_data_nojuly$District, 
                          k = minlag:maxlag)
lag_hh = tsModel::Lag(monthly_dengue_data_nojuly$Average_Household, 
                      group = monthly_dengue_data_nojuly$District, 
                      k = minlag:maxlag)

lag_precip_cut <- lag_precip[monthly_dengue_data_nojuly$year > 2013,]
lag_tmax_cut <- lag_tmax[monthly_dengue_data_nojuly$year > 2013,]
lag_tmin_cut <- lag_tmin[monthly_dengue_data_nojuly$year > 2013,]
lag_hum_cut <- lag_hum[monthly_dengue_data_nojuly$year > 2013,]
lag_tcases_cut <- lag_tcases[monthly_dengue_data_nojuly$year > 2013,]
lag_popden_cut <- lag_popden[monthly_dengue_data_nojuly$year > 2013,]
lag_hh_cut <- lag_hh[monthly_dengue_data_nojuly$year > 2013,]

lagknot = equalknots(3, 12, 2) 

monthly_dengue_data_nojuly <- monthly_dengue_data_nojuly[monthly_dengue_data_nojuly$year > 2013,]
monthly_dengue_data <- monthly_dengue_data[monthly_dengue_data$year > 2013,]

# construct cross basis matrices

var <- lag_precip_cut
basis_precip_cut <- crossbasis(var,
                               argvar = list(fun = 'ns', 
                                             knots = equalknots(monthly_dengue_data_nojuly$total_precip, 2)),
                               arglag = list(fun = 'ns', 
                                             knots = (maxlag-minlag)/2))

var <- lag_tmax_cut
basis_tmax_cut <- crossbasis(var,
                             argvar = list(fun = 'ns', 
                                           knots = equalknots(monthly_dengue_data_nojuly$avg_tempmax, 2)),
                             arglag = list(fun = 'ns', 
                                           knots = (maxlag-minlag)/2))

var <- lag_tmin_cut
basis_tmin_cut <- crossbasis(var,
                             argvar = list(fun = 'ns', 
                                           knots = equalknots(monthly_dengue_data_nojuly$avg_tempmin, 2)),
                             arglag = list(fun = 'ns', 
                                           knots = (maxlag-minlag)/2))

var <- lag_hum_cut
basis_hum_cut <- crossbasis(var,
                            argvar = list(fun = 'ns', 
                                          knots = equalknots(monthly_dengue_data_nojuly$avg_humidity, 2)),
                            arglag = list(fun = 'ns', 
                                          knots = (maxlag-minlag)/2))

var <- lag_tcases_cut
basis_tcases_cut <- crossbasis(var,
                               argvar = list(fun = 'ns', 
                                             knots = equalknots(monthly_dengue_data_nojuly$total_cases, 2)),
                               arglag = list(fun = 'ns', 
                                             knots = (maxlag-minlag)/2))

var <- lag_popden_cut
basis_popden_cut <- crossbasis(var,
                               argvar = list(fun = 'ns', 
                                             knots = equalknots(monthly_dengue_data_nojuly$PopDen, 2)),
                               arglag = list(fun = 'ns', 
                                             knots = (maxlag-minlag)/2))

var <- lag_hh_cut
basis_hh_cut <- crossbasis(var,
                           argvar = list(fun = 'ns', 
                                         knots = equalknots(monthly_dengue_data_nojuly$Average_Household, 2)),
                           arglag = list(fun = 'ns', 
                                         knots = (maxlag-minlag)/2))
View(basis_precip_cut)

colnames(basis_precip_cut) = paste0('b_precip', colnames(basis_precip_cut))
colnames(basis_tmax_cut) = paste0('b_tmax', colnames(basis_tmax_cut))
colnames(basis_tmin_cut) = paste0('b_tmin', colnames(basis_tmin_cut))
colnames(basis_hum_cut) = paste0('b_hum', colnames(basis_hum_cut))
colnames(basis_tcases_cut) = paste0('b_tcases', colnames(basis_tcases_cut))
colnames(basis_popden_cut) = paste0('b_popden', colnames(basis_popden_cut))
colnames(basis_hh_cut) = paste0('b_hh', colnames(basis_hh_cut))


############################################################################### model specification

# model with cross basis (MODEL A)
formula_cross <- total_cases ~ basis_tcases_cut + 
  basis_precip_cut +
  basis_tmax_cut +
  basis_tmin_cut +
  basis_hum_cut +
  basis_hh_cut +
  basis_popden_cut +
  f(year, month, model = "rw2") + 
  f(District, model = "besag", graph = as(as.matrix(st_touches(districts)), "CsparseMatrix"))

# reduced model (MODEL B)
formula1 <- total_cases ~ total_cases_lag_1 + lag_precip_2 + lag_tmin_2 + lag_hum_2 +
  f(year, month, model = "rw2") +
  f(District, model = "besag", graph = as(as.matrix(st_touches(districts)), "CsparseMatrix"))

# model without intervention (MODEL C)
formula2 <- total_cases ~ total_cases_lag_1 + lag_precip_2 + lag_tmax_2 + lag_tmin_2 + lag_hum_2 + hh_lag_1+
  popden_lag_1  +
  f(year, month, model = "rw2") + 
  f(District, model = "besag", graph = as(as.matrix(st_touches(districts)), "CsparseMatrix"))


############################################################################### model fitting

# MODEL A
## fit the spatiotemporal model using INLA with a negative binomial distribution
spatiotemporal_model_cross <- inla(formula_cross,
                                   data = monthly_dengue_data_nojuly,
                                   family = 'nbinomial',
                                   control.predictor = list(link = 1, compute = TRUE),
                                   control.inla = list(strategy = 'adaptive'),
                                   control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE))

# MODEL B
spatiotemporal_model1 <- inla(formula1,
                              data = monthly_dengue_data_nojuly,
                              family = "nbinomial",
                              control.predictor = list(compute = TRUE),
                              control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config=TRUE))

# MODEL C
spatiotemporal_model2 <- inla(formula2, 
                              data = monthly_dengue_data_nojuly, 
                              family = "nbinomial", 
                              control.predictor = list(compute = TRUE), 
                              control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config=TRUE))


############################################################################### model metrics

# function that prints out model's metrics (DIC, WAIC, Log Score)
model_metrics <- function(model) {
  dic <- model$dic$dic
  waic <- model$waic$waic
  cpo_values <- model$cpo$cpo
  log_score <- sum(log(cpo_values))
  print(paste("DIC:", dic))
  print(paste("WAIC:", waic))
  print(paste("Logarithmic Score of CPO:", log_score))
}

# MODEL A
model_metrics(spatiotemporal_model_cross)

# [1] "DIC: 6880.10981991521"
# [1] "WAIC: 6891.9794410332"

# MODEL B
model_metrics(spatiotemporal_model1)

# [1] "DIC: 7328.0894209513"
# [1] "WAIC: 7330.25028453821"

# MODEL C
model_metrics(spatiotemporal_model2)

# [1] "DIC: 7323.13315898655"
# [1] "WAIC: 7325.71063417704"


############################################################################### exceedance probability visualizations

exceedance_viz <- function(model, threshold) {
  
  # find fitted values obtained from model
  fitted_values <- model$summary.fitted.values$mean
  
  # set threshold 
  threshold = threshold
  
  # get estimated overdispersion parameter
  size <- model$summary.hyperpar$mean[1]
  
  # calculate exceedance probability
  exceedance_prob <- 1 - pnbinom(threshold, size = size, mu = fitted_values)
  
  # add exceedance probabilities directly to dataset
  monthly_dengue_data_nojuly$exceedance_prob <- exceedance_prob  
  
  # create a factor for months to ensure correct ordering in plots
  monthly_dengue_data_nojuly$month_f <- factor(monthly_dengue_data_nojuly$month, levels = 1:12, labels = month.abb)
  
  # convert data to sf object
  dengue_cases_vis <- monthly_dengue_data_nojuly[, c("year", "month", "District", "exceedance_prob","longitude","latitude")]
  dengue_cases_current <- monthly_dengue_data_nojuly[, c("year", "month", "District", "total_cases","longitude","latitude")]
  
  # merge the dengue cases data with the districts spatial data
  merged_data <- st_join(districts, dengue_cases_vis)
  merged_current_data <- st_join(districts, dengue_cases_current)
  
  # filter data for year and month
  merged_data_2024 <- merged_data %>%
    filter(year == 2024, month == 7)
  
  exceedance_vis <- ggplot(merged_data_2024) +
    geom_sf(aes(fill = exceedance_prob), color = "black") +
    scale_fill_gradient(name = "Exceedance Probability", 
                        low = "lightblue", high = "darkblue",
                        guide = guide_colorbar(title.position = "top", 
                                               title.hjust = 0.5,
                                               label.position = "bottom",
                                               label.hjust = 0.5,
                                               direction = "horizontal")) +
    labs(title = "Exceedance Probability ",
         subtitle = "of Dengue Cases on July 2024",
         caption = "(Threshold = 78) ") +
    theme_minimal() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 12,face = "bold"),
          plot.caption = element_text(size = 10, hjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 5),
          axis.text = element_blank(),
          axis.title = element_blank())
  print(exceedance_vis)
  
}

# MODEL A
exceedance_viz(spatiotemporal_model_cross, 78)

# MODEL B
exceedance_viz(spatiotemporal_model1, 78)

# MODEL C
exceedance_viz(spatiotemporal_model2, 78)


############################################################################### prediction graphs

prediction_viz <- function(model) {
  
  # extract summary fitted values
  summary_fitted_values <- model$summary.fitted.values
  
  # extract posterior mean, standard deviation, and credible intervals (2.5%, 97.5%)
  predicted_mean <- summary_fitted_values$mean
  predicted_sd <- summary_fitted_values$sd
  predicted_lower <- summary_fitted_values$`0.025quant`
  predicted_upper <- summary_fitted_values$`0.975quant`
  
  monthly_dengue_data_nojuly <- monthly_dengue_data_nojuly %>%
    mutate(predicted_mean = predicted_mean,
           predicted_sd = predicted_sd,
           predicted_lower = predicted_lower,
           predicted_upper = predicted_upper)

  # filter dataset by specific month and year
  filtered_data <- monthly_dengue_data_nojuly %>%
    filter(year == 2024, month == 7)
  
  # plot predictions
  #ggplot(filtered_data, aes(x = District, y = predicted_mean_c)) +
  #  geom_bar(stat = "identity", fill = "blue") +
  #  labs(title = "Predicted Dengue Cases for July 2024",
  #       x = "District",
  #       y = "Predicted Cases") +
  #  theme_minimal()
  
  # plot predictions with points and error bars
  ggplot(filtered_data, aes(x = District, y = predicted_mean)) +
    geom_point(size = 3, color = "blue") +  # Plot the points for the predicted means
    geom_errorbar(aes(ymin = predicted_lower, ymax = predicted_upper), width = 0.2) +  # Add error bars
    labs(title = "Predicted Dengue Cases for July 2024",
         x = "District",
         y = "Predicted Cases") +
    ylim(00, 600) +
    theme_minimal() 
  
}

# MODEL A
prediction_viz(spatiotemporal_model_cross)

# MODEL B
prediction_viz(spatiotemporal_model1)

# MODEL C
prediction_viz(spatiotemporal_model2)


############################################################################### comparison graphs

comparison_viz <- function(model) {
  
  # extract summary fitted values
  summary_fitted_values <- model$summary.fitted.values
  
  # extract posterior mean, standard deviation, and credible intervals (2.5%, 97.5%)
  predicted_mean <- summary_fitted_values$mean
  predicted_sd <- summary_fitted_values$sd
  predicted_lower <- summary_fitted_values$`0.025quant`
  predicted_upper <- summary_fitted_values$`0.975quant`
  
  monthly_dengue_data <- monthly_dengue_data %>%
    mutate(predicted_mean = predicted_mean,
           predicted_sd = predicted_sd,
           predicted_lower = predicted_lower,
           predicted_upper = predicted_upper)
  
  filtered_data <- monthly_dengue_data %>%
    filter(year == 2024, month == 7)
  
  filtered_data_long <- filtered_data %>%
    select(District, predicted_mean, total_cases) %>%
    pivot_longer(cols = c(predicted_mean, total_cases), names_to = "Type", values_to = "Cases")
  
  ggplot(filtered_data_long, aes(x = District, y = Cases, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Predicted vs Actual Dengue Cases for July 2024",
         x = "District",
         y = "Number of Cases",
         fill = "Case Type") +
    scale_fill_manual(values = c("predicted_mean" = "blue", "total_cases" = "red")) +
    theme_minimal() +
    ylim(0,250)
}

comparison_line <- function(model) {
  
  # extract summary fitted values
  summary_fitted_values <- model$summary.fitted.values
  
  # extract posterior mean, standard deviation, and credible intervals (2.5%, 97.5%)
  predicted_mean <- summary_fitted_values$mean
  predicted_sd <- summary_fitted_values$sd
  predicted_lower <- summary_fitted_values$`0.025quant`
  predicted_upper <- summary_fitted_values$`0.975quant`
  
  monthly_dengue_data <- monthly_dengue_data %>%
    mutate(predicted_mean = predicted_mean,
           predicted_sd = predicted_sd,
           predicted_lower = predicted_lower,
           predicted_upper = predicted_upper)
  
  filtered_data_line <- monthly_dengue_data %>%
    mutate(Date = as.Date(paste(year, month, "01", sep = "-"), format = "%Y-%m-%d"))
  
  ggplot(filtered_data_line, aes(x = Date)) +
    geom_ribbon(aes(ymin = predicted_lower, ymax = predicted_upper), fill = "gray80", alpha = 0.5) +
    geom_line(aes(y = total_cases, color = "Cases")) +
    geom_line(aes(y = predicted_mean, color = "Predicted Mean")) +
    facet_wrap(~ District, scales = "free_y") +
    labs(title = "Dengue Cases and Predicted Mean over Time ",
         x = "Time",
         y = "Number of Cases",
         color = "Legend" ) + 
    theme_minimal()
  
}

# MODEL A
comparison_viz(spatiotemporal_model_cross)

# MODEL B
comparison_viz(spatiotemporal_model1)

# MODEL C
comparison_viz(spatiotemporal_model2)



# MODEL A
comparison_line(spatiotemporal_model_cross)

# MODEL B
comparison_line(spatiotemporal_model1)

# MODEL C
comparison_line(spatiotemporal_model2)
