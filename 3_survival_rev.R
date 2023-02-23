# LOAD REQUIRED PACKAGES ####
library(lubridate)
library(MASS)
library(bbmle)
library(lsmeans)
library(ggplot2)
library(patchwork)

# LOAD DATA ####
# ensure PIT tag numbers do not present in scientific notation
options(scipen=999)

# load comprehensive mark data set (all fish tagged at RST with predicted length at marine entry)
### individual fish are identified by PIT tag number
dat_rst <-read.csv("dat_rst_mar.csv") 

# load comprehensive mark-recapture data set (only tagged fish that were recaptured in the estuary)
### individual fish are identified by PIT tag number
### columns preceded by "tag_" indicate initial mark data
### columns preceded by "rc_" indicate recapture data
dat_rc <- read.csv("dat_rc.csv")

# load comprehensive marine survival data set (all fish tagged at RST with adult survival data from RFID array detections)
### 0 = died, 1 = redetected as adult 
dat_surv <- read.csv("dat_surv.csv")

# CHECK: SIZE-SELECTIVE MORTALITY BETWEEN FRESHWATER EXIT AND MARINE ENTRY ####
### annual logistic regressions to test whether length at tagging significantly different for tagged at the RST vs those recaptured in Reach 1 
### no significant difference (slope not significantly different from 0; p > 0.05) inidicates a lack of size-selective mortality in the estuary

### create annual RST tagging data frames
dat_rst_17 <- subset(dat_rst, tag_year==2017)
dat_rst_18 <- subset(dat_rst, tag_year==2018)
dat_rst_19 <- subset(dat_rst, tag_year==2019)

### create annual Reach 1 recapture data frames
dat_rc <- subset(dat_rc, tag_loc == "smolt trap" & rc_loc_1 == 1) # fish tagged at RST and recaptured in Reach 1 (estuary mouth / marine entry)
dat_rc$tag_year <- as.factor(dat_rc$tag_year)
dat_rc$tag_year <- droplevels(dat_rc$tag_year)

dat_rc_17 <- subset(dat_rc, tag_year==2017) # Reach 1 recaptures in 2017
dat_rc_18 <- subset(dat_rc, tag_year==2018) # Reach 1 recaptures in 2018
dat_rc_19 <- subset(dat_rc, tag_year==2019) # Reach 1 recaptures in 2019

### 2017 logistic regression: no significant effect of length; p = 0.16770
dat_rst_17$rc <- 0
dat_rst_17$rc[match(dat_rc_17$PIT_number, dat_rst_17$PIT_number)] <- 1
bin_17 <- glm(rc ~ tag_length, data = dat_rst_17, family = binomial("logit"));summary(bin_17) 

# 2018 logistic regression: no significant effect of length; p = 0.53352
dat_rst_18$rc <- 0
dat_rst_18$rc[match(dat_rc_18$PIT_number, dat_rst_18$PIT_number)] <- 1
bin_18 <- glm(rc ~ tag_length, data = dat_rst_18, family = binomial("logit"));summary(bin_18) 

### 2019 no significant effect of length; p = 0.973
dat_rst_19$rc <- 0
dat_rst_19$rc[match(dat_rc_19$PIT_number, dat_rst_19$PIT_number)] <- 1
bin_19 <- glm(rc ~ tag_length, data = dat_rst_19, family = binomial("logit"));summary(bin_19) 

# PREPARE MARINE SURVIVAL DATA FOR SURVIVAL MODELS ####
# keep only columns relevant for analyses
dat_surv <- dat_surv[,c("PIT_number", "tag_year", "doy", "length", "surv", "recap_year")]
colnames(dat_surv) <- c("PIT_number", "tag_year", "tag_doy", "length_fw", "surv", "rc_year")

# check number of tagged spawners returning each year
table(dat_surv$tag_year, dat_surv$rc_year)

# append marine entry sizes to marine survival data
### create columns for predicted length at marine entry (length_mar) and standard deviation of growth estimate (sd_growth)
dat_surv$length_mar <- NA
dat_surv$sd_growth <- NA

# append length at marine entry and growth SD values to survival data frame
### for each individual in the survival data set
for(i in 1:nrow(dat_surv)) {
  ### identify tag year 
  tag_year <- dat_surv$tag_year[i]
  
  ### identify length at freshwater exit (tag length)
  tag_length <- dat_surv$length_fw[i]
  
  ### identify fish in RST data set that matches tag year and tag length
  id <- which(dat_rst$tag_year == tag_year & dat_rst$tag_length == tag_length)
  
  ### for that fish, assign predicted length at marine entry
  dat_surv$length_mar[i] <- unique(dat_rst[id,]$pred_length)
  
  ### for that fish, assign growth variance
  dat_surv$var_growth[i] <- unique(dat_rst[id,]$pred_growth_var)
}

# restrict survival data to fish tagged in 2017 or 2018
### these are smolt years with complete adult return data available
dat_surv <- dat_surv[dat_surv$tag_year == 2017 | dat_surv$tag_year == 2018,]

# create two data frames for marine survival models
### both data frames include all fish tagged at the RST in 2017 and 2018
### one data frame is for observed survival model (estuary-present model) -- length is observed length at freshwater exit; estuary growth is implicit in survival data because fish transit the estuary en route to the ocean
### one data frame is for predicted survival model (estuary-absent model) -- length is modeled length at marine entry; no estuary growth between freshwater exit and marine entry b/c these points are coincident in space and time
surv_obs <- data.frame("tag_year" = as.factor(dat_surv$tag_year),
                      "length" = dat_surv$length_fw,
                      "surv" = dat_surv$surv,
                      "estuary" = 1,
                      "var_growth" = NA,
                      "weight" = NA,
                      "weight_sc" = 1 # all fish have equal weight in the model
                      )

surv_pred <- data.frame("tag_year" = as.factor(dat_surv$tag_year),
                        "length" = dat_surv$length_mar, 
                        "surv" = dat_surv$surv,
                        "estuary" = 0,
                        "var_growth" = dat_surv$var_growth,# variance of predicted days + variance of predicted daily growth - covariance between predicted days and predicted daily growth (per reviewer suggestion)
                        "weight" = NA, 
                        "weight_sc" = NA
                        )

# for estuary-absent model, propagate error from residence and growth model predictions to weight individual fish
### each individual's weight is the inverse of the variance associated with its predicted size at marine entry -- this propagates error associated with both estuary residence and growth models per revierwer comments, such that fish with more variance about their predicted length receive less weight in model fitting
### weights are scaled to 1 for consistency with observed survival data
surv_pred$weight <- 1/(surv_pred$var_growth^2) 
surv_pred$weight_sc <- surv_pred$weight/max(surv_pred$weight) 

# combine observed and predicted data sets into one data frame
dat <- data.frame(rbind(surv_obs, surv_pred))
dat$surv <- as.integer(dat$surv)
dat$estuary <- as.factor(dat$estuary)
dat$tag_year <- as.factor(dat$tag_year)

# RUN CANDIDATE MARINE SURVIVAL MODELS ####
# RESPONSE VARAIABLE: observed survival (surv)
# CANDIDATE PREDICTORS: length (length), tag year (tag_year), and estuary (1 = estuary present, 0 = estuary absent)
# MODEL STRUCTURE: binomial (0/1 survival data)

### run candidate models
mod_null <- glm(surv ~ 1, family = binomial("logit"), data = dat, weights = weight_sc);summary(mod_null) # null model - constant survival rate

mod_l <- glm(surv ~ length, family = binomial("logit"), data = dat, weights = weight_sc);summary(mod_l) # effect of length

mod_l_y <- glm(surv ~ length + tag_year, family = binomial("logit"), data = dat, weights = weight_sc);summary(mod_l_y) # effect of length + tag year

mod_l_y_est <- glm(surv ~ length + tag_year + estuary, family = binomial("logit"), data = dat, weights = weight_sc);summary(mod_l_y_est) # effect of length + tag year + estuary

mod_l_y_est_int <- glm(surv ~ length + tag_year*estuary, family = binomial("logit"), data = dat, weights = weight_sc);summary(mod_l_y_est_int)

### use AICc to calculate best-fit model for binomial distribution
# determine best-fit annual growth rate model
### aggregate candidate models into list
mod_list <- list("mod_null" = mod_null, # null model - constant growth rate
                 "mod_l" = mod_l, # effect of tag length
                 "mod_l_y" = mod_l_y, # effect of tag length + tag year
                 "mod_l_y_est" = mod_l_y_est, # effect of tag length + tag year + categorical estuary
                 "mod_l_y_est_int" = mod_l_y_est_int # effect of tag length + tag_year + estuary:year
)

# TABLE 3: MARINE SURVIVAL MODEL SELECTION ####
### compare models using AICc
AICtab(mod_list, delta = TRUE, weights = TRUE, base = TRUE)

### best-fit marine survival model: surv ~ length + tag year + estuary
summary(mod_l_y_est)


# PAIRWISE ANNUAL SURVIVAL COMPARISONS ####
#  tukey test to determine whether annual survival rates are significantly different from one another
### 2017 differs from 2018
mod_surv_lsm <- lsmeans(mod_l_y_est, ~tag_year)
contrast(mod_surv_lsm, "tukey")

# PREDICT MARINE SURVIVAL PROBABILITY FOR EACH FISH TAGGED AT THE RST #### 
# standardized lenght inputs for model predictions
length_vals <- seq(min(dat_surv$length_fw[dat_surv$surv == 1]), max(dat_surv$length_fw[dat_surv$surv == 1]), 0.1)

# categorical estuary and year inputs for model predictions
est_0 <- as.factor(rep(0, length(length_vals)))
est_1 <- as.factor(rep(1, length(length_vals)))
year_17 <- as.factor(rep(2017, length(length_vals)))
year_18 <- as.factor(rep(2018, length(length_vals)))

# predict survival probabilities; create data frames with CIs
### 2017 estuary-present
surv_obs_17 <- predict(mod_l_y_est, list(length = length_vals, estuary = est_1, tag_year = year_17), se.fit = T, type = "link")
surv_obs_17 <- as.data.frame(surv_obs_17)
surv_obs_17$lower <- surv_obs_17$fit - (1.96*surv_obs_17$se.fit)
surv_obs_17$upper <- surv_obs_17$fit + (1.96*surv_obs_17$se.fit)

surv_obs_17$fit_odds <- exp(surv_obs_17$fit)
surv_obs_17$fit_prob <- surv_obs_17$fit_odds/(1+surv_obs_17$fit_odds)

surv_obs_17$lower_odds <- exp(surv_obs_17$lower)
surv_obs_17$lower_prob <- surv_obs_17$lower_odds/(1+surv_obs_17$lower_odds)

surv_obs_17$upper_odds <- exp(surv_obs_17$upper)
surv_obs_17$upper_prob <- surv_obs_17$upper_odds/(1+surv_obs_17$upper_odds)

### 2018 estuary-present
surv_obs_18 <- predict(mod_l_y_est, list(length = length_vals, estuary = est_1, tag_year = year_18), se.fit = T, type = "link")
surv_obs_18 <- as.data.frame(surv_obs_18)
surv_obs_18$lower <- surv_obs_18$fit - (1.96*surv_obs_18$se.fit) 
surv_obs_18$upper <- surv_obs_18$fit + (1.96*surv_obs_18$se.fit) 

surv_obs_18$fit_odds <- exp(surv_obs_18$fit)
surv_obs_18$fit_prob <- surv_obs_18$fit_odds/(1+surv_obs_18$fit_odds)

surv_obs_18$lower_odds <- exp(surv_obs_18$lower)
surv_obs_18$lower_prob <- surv_obs_18$lower_odds/(1+surv_obs_18$lower_odds)

surv_obs_18$upper_odds <- exp(surv_obs_18$upper)
surv_obs_18$upper_prob <- surv_obs_18$upper_odds/(1+surv_obs_18$upper_odds)


### 2017 estuary-absent
surv_pred_17 <- predict(mod_l_y_est, list(length = length_vals, estuary = est_0, tag_year = year_17), se.fit = T, type = "link")
surv_pred_17 <- as.data.frame(surv_pred_17)
surv_pred_17$lower <- surv_pred_17$fit - (1.96*surv_pred_17$se.fit) 
surv_pred_17$upper <- surv_pred_17$fit + (1.96*surv_pred_17$se.fit) 

surv_pred_17$fit_odds <- exp(surv_pred_17$fit)
surv_pred_17$fit_prob <- surv_pred_17$fit_odds/(1+surv_pred_17$fit_odds)

surv_pred_17$lower_odds <- exp(surv_pred_17$lower)
surv_pred_17$lower_prob <- surv_pred_17$lower_odds/(1+surv_pred_17$lower_odds)

surv_pred_17$upper_odds <- exp(surv_pred_17$upper)
surv_pred_17$upper_prob <- surv_pred_17$upper_odds/(1+surv_pred_17$upper_odds)


### 2018 estuary-absent
surv_pred_18 <- predict(mod_l_y_est, list(length = length_vals, estuary = est_0, tag_year = year_18), se.fit = T, type = "link")
surv_pred_18 <- as.data.frame(surv_pred_18)
surv_pred_18$lower <- surv_pred_18$fit - (1.96*surv_pred_18$se.fit) 
surv_pred_18$upper <- surv_pred_18$fit + (1.96*surv_pred_18$se.fit) 

surv_pred_18$fit_odds <- exp(surv_pred_18$fit)
surv_pred_18$fit_prob <- surv_pred_18$fit_odds/(1+surv_pred_18$fit_odds)

surv_pred_18$lower_odds <- exp(surv_pred_18$lower)
surv_pred_18$lower_prob <- surv_pred_18$lower_odds/(1+surv_pred_18$lower_odds)

surv_pred_18$upper_odds <- exp(surv_pred_18$upper)
surv_pred_18$upper_prob <- surv_pred_18$upper_odds/(1+surv_pred_18$upper_odds)


# SAVE DATA, MODELS, AND PREDICTIONS ####
saveRDS(dat, "dat_surv_obs_pred.rds")
saveRDS(mod_l_y_est, "mod_surv.rds")
saveRDS(surv_obs_17, "surv_obs_17.rds")
saveRDS(surv_obs_18, "surv_obs_18.rds")
saveRDS(surv_pred_17, "surv_pred_17.rds")
saveRDS(surv_pred_18, "surv_pred_18.rds")

# TABLE 4: MARINE SURVIVAL MODEL EFFECT SIZES, CONFIDENCE INTERVALS, P-VALUES ####
# run model without intercept (for interpretability)
mod_l_y_est <- glm(surv ~ 0 + length + tag_year + estuary, family = binomial("logit"), data = dat, weights = weight_sc)

# determine effect sizes, p-values
summary(mod_l_y_est)

# determine confidence intervals
round(confint(mod_l_y_est), 3)

# FIGURE 6: MARINE SURVIVAL ####
# create survival data frames for easy plotting
### estuary-present model
p_surv_obs <- rbind(surv_obs_17, surv_obs_18)
p_surv_obs$length <- rep(length_vals,2)
p_surv_obs$year <- c(as.character(year_17), as.character(year_18))

### estuary-absent model
p_surv_pred <- rbind(surv_pred_17, surv_pred_18)
p_surv_pred$length <- rep(length_vals,2)
p_surv_pred$year <- c(as.character(year_17), as.character(year_18))

# establish plot inputs
### points
pt_size <- 1.5
pt_co <- "grey50"

### lines
line_size <- 1
obs_co <- "#023a45"
pred_co <- "#2baa8e"

### axes
ax_title <- 10
ax_text <- 8

### quantile points
quantile(dat_rst$tag_length, c(0.1, 0.5, 0.9))
q0.1 <- 78
q0.5 <- 93
q0.9 <- 112

# determine predicted decrease in survival for median fish w/o estuary rearing
obs0.5_17 <- p_surv_obs[p_surv_obs$length==q0.5 & p_surv_obs$year==2017,]
pred0.5_17 <- p_surv_pred[p_surv_pred$length==q0.5 & p_surv_pred$year==2017,]
pred0.5_17$length_pred <- p_surv_obs$length[match(round(pred0.5_17$fit_prob, 4), round(p_surv_obs$fit_prob[p_surv_obs$year==2017],4))]

obs0.5_18 <- p_surv_obs[p_surv_obs$length==q0.5 & p_surv_obs$year==2018,]
pred0.5_18 <- p_surv_pred[p_surv_pred$length==q0.5 & p_surv_pred$year==2018,]
pred0.5_18$length_pred <- p_surv_obs$length[match(round(pred0.5_18$fit_prob, 4), round(p_surv_obs$fit_prob[p_surv_obs$year==2018],4))]

# predicted increase in survival for median fish w/o estuary rearing
p_surv_obs$p[p_surv_obs$year==2017] <- (p_surv_obs$fit_prob[p_surv_obs$year==2017]/p_surv_pred$fit_prob[p_surv_pred$year==2017])-1
p_surv_obs$p_lower[p_surv_obs$year==2017] <- (p_surv_obs$lower_prob[p_surv_obs$year==2017]/p_surv_pred$lower_prob[p_surv_pred$year==2017])-1
p_surv_obs$p_upper[p_surv_obs$year==2017] <- (p_surv_obs$upper_prob[p_surv_obs$year==2017]/p_surv_pred$upper_prob[p_surv_pred$year==2017])-1

p_surv_obs$p[p_surv_obs$year==2018] <- (p_surv_obs$fit_prob[p_surv_obs$year==2018]/p_surv_pred$fit_prob[p_surv_pred$year==2018])-1
p_surv_obs$p_lower[p_surv_obs$year==2018] <- (p_surv_obs$lower_prob[p_surv_obs$year==2018]/p_surv_pred$lower_prob[p_surv_pred$year==2018])-1
p_surv_obs$p_upper[p_surv_obs$year==2018] <- (p_surv_obs$upper_prob[p_surv_obs$year==2018]/p_surv_pred$upper_prob[p_surv_pred$year==2018])-1

p_surv_obs[p_surv_obs$length==q0.1 & p_surv_obs$year==2017,]
p_surv_obs[p_surv_obs$length==q0.5 & p_surv_obs$year==2017,]
p_surv_obs[p_surv_obs$length==q0.9 & p_surv_obs$year==2017,]

p_surv_pred[p_surv_pred$length==q0.1 & p_surv_pred$year==2017,]
p_surv_pred[p_surv_pred$length==q0.5 & p_surv_pred$year==2017,]
p_surv_pred[p_surv_pred$length==q0.9 & p_surv_pred$year==2017,]

p_surv_obs$fit_prob[p_surv_obs$length==q0.1 & p_surv_obs$year==2017]/p_surv_pred$fit_prob[p_surv_pred$length==q0.1 & p_surv_pred$year==2017]
p_surv_obs$fit_prob[p_surv_obs$length==q0.5 & p_surv_obs$year==2017]/p_surv_pred$fit_prob[p_surv_pred$length==q0.5 & p_surv_pred$year==2017]
p_surv_obs$fit_prob[p_surv_obs$length==q0.9 & p_surv_obs$year==2017]/p_surv_pred$fit_prob[p_surv_pred$length==q0.9 & p_surv_pred$year==2017]

p_surv_obs[p_surv_obs$length==q0.1 & p_surv_obs$year==2018,]
p_surv_obs[p_surv_obs$length==q0.5 & p_surv_obs$year==2018,]
p_surv_obs[p_surv_obs$length==q0.9 & p_surv_obs$year==2018,]

p_surv_pred[p_surv_pred$length==q0.1 & p_surv_pred$year==2018,]
p_surv_pred[p_surv_pred$length==q0.5 & p_surv_pred$year==2018,]
p_surv_pred[p_surv_pred$length==q0.9 & p_surv_pred$year==2018,]

p_surv_obs$fit_prob[p_surv_obs$length==q0.1 & p_surv_obs$year==2018]/p_surv_pred$fit_prob[p_surv_pred$length==q0.1 & p_surv_pred$year==2018]
p_surv_obs$fit_prob[p_surv_obs$length==q0.5 & p_surv_obs$year==2018]/p_surv_pred$fit_prob[p_surv_pred$length==q0.5 & p_surv_pred$year==2018]
p_surv_obs$fit_prob[p_surv_obs$length==q0.9 & p_surv_obs$year==2018]/p_surv_pred$fit_prob[p_surv_pred$length==q0.9 & p_surv_pred$year==2018]

# draw figure 6: marine survival plot
fig_6a <- ggplot(p_surv_obs) + 
  facet_grid(year ~ ., switch = "y") + 
  coord_cartesian(xlim = c(60, 142), ylim = c(0, 0.18)) +
  scale_x_continuous(name = "Fork length at at freshwater exit (mm)", breaks = seq(60, 140, 20)) +
  scale_y_continuous(
  
    # Features of the primary y-axis
    name = "Probability of survival", breaks = seq(0, 0.16, 0.04),
    
    # Features of the secondary y-axis
    sec.axis = sec_axis( trans=~.*5.55, name="Absolute survival", breaks = c(0,1))
  ) +
  
  # estuary-absent model predictions
  geom_line(data = p_surv_pred, aes(x = length, y = fit_prob, color = "white"), lwd=line_size) +

  # estuary-present model predictions
  geom_line(data = p_surv_obs, aes(x = length, y = fit_prob, color = obs_co), lwd=line_size) +
  geom_ribbon(data = p_surv_obs, aes(ymin = lower_prob, ymax = upper_prob, x = length), fill = adjustcolor(obs_co, 0.4)) +
  
  # add points for median fish - observed
  geom_point(data = p_surv_obs[p_surv_obs$length==q0.5,], aes(x = length, y = fit_prob), size = 2.5, color = obs_co, fill = obs_co) +

  # add points for median fish - predicted
  geom_point(data = pred0.5_17, aes(x = length_pred, y = fit_prob), size = 2.5, color = pred_co, fill = pred_co) +
  geom_point(data = pred0.5_18, aes(x = length_pred, y = fit_prob), size = 2.5, color = pred_co, fill = pred_co) +
  
  # add arrows showing decrease in survival probability
  geom_curve(data = obs0.5_17, aes(x = length, y = fit_prob + 0.013, xend = pred0.5_17$length_pred-2, yend = pred0.5_17$fit_prob + 0.010), curvature = -0.031, ncp = 10, arrow=arrow(angle = 15, ends = "both", length = unit(0.1, "inches"), type = "closed"), color = obs_co, size = line_size*0.75, lty = 1)+
  geom_curve(data = obs0.5_18, aes(x = length, y = fit_prob + 0.015, xend = pred0.5_18$length_pred-2, yend = pred0.5_18$fit_prob + 0.012), curvature = -0.031, ncp = 10, arrow=arrow(angle = 15, ends = "both", length = unit(0.1, "inches"), type = "closed"), color = obs_co, size = line_size*0.75, lty = 1)+
  
  # add secondary axis points indicating absolute individual survival (0/1) 
  geom_point(data = dat_surv[dat_surv$surv==0,], aes(x = length_fw, y = 0, alpha = 0.6), position = position_jitter(w = 1, h = 0), col = pt_co, shape = 22, size = 2.5, show.legend = F) +
  geom_point(data = dat_surv[dat_surv$surv==1,], aes(x = length_fw, y = 0.18, alpha = 0.6), position = position_jitter(w = 1, h = 0), col = pt_co, shape = 22, size = 2.5, show.legend = F) +
  
  scale_color_identity(name = "Marine survival model",
                       breaks = c(obs_co, "white"),
                       labels = c("Estuary-present model (obs)",
                                  "Estuary-absent model (pred)"),
                       guide = "legend") +

  
  # map theme
  theme_classic() +
  theme(axis.title = element_text(size = ax_title, face = "bold"), 
        axis.text = element_text(size = ax_text, color = "black"),
        axis.line = element_line(color = NA),
        strip.background = element_rect(fill = "grey70"),
        strip.text = element_text(size = ax_title, face = "bold"),
        panel.border = element_rect(fill = NA, size = 1, color = "black"),
        legend.title = element_text(size = ax_text, color = "black"),
        legend.text = element_text(size = ax_text-0.65, color = "black"),
        legend.position = c(0.45, 0.87))

fig_6a

# subset data for figure 6b
dat_p <- subset(p_surv_obs, length %in% c(q0.1, q0.5, q0.9))

fig_6b <- ggplot(dat_p) +
  facet_grid(year ~ .) + 
  coord_cartesian(ylim = c(0,80)) +
  scale_x_discrete(labels = c("Q0.1", "Q0.5","Q0.9")) +
  scale_y_continuous(name = "Increase in probability of survival (%)", position = "right") +
  geom_linerange(aes(x = factor(length), ymin = p_lower*100, ymax = p_upper*100), position = position_dodge(0.5), size = line_size, color = obs_co) +
  geom_point(aes(x = factor(length), y = p*100, shape = factor(length)), position= position_dodge(width = 0.5), size = 2.5, color = obs_co, fill = obs_co) +
  scale_shape_manual(values = c(23, 21, 22), labels = c("Q0.1", "Q0.5","Q0.9")) +
  labs(x = "Size quantile") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = ax_title, face = "bold"), 
        axis.text = element_text(size = ax_text, color = "black"),
        axis.line = element_line(color = NA),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(fill = NA, size = 1, color = "black"))

fig_6b

# save marine survival plot
(fig_6a + fig_6b) + plot_layout(widths = c(2,1))
ggsave("fig_6_rev.pdf", dpi = 300, height = 6, width = 6, units = "in")


# REVEIEWER COMMENT RESPONSE - SIZE AT MARINE ENTRY VALIDATION PLOT ####
### plot observed vs. predicted size at marine entry (reafch 2)
### all fish observed at the RST
dat_rst_val <- dat_rst[,c("tag_length", "pred_length", "tag_year")]
colnames(dat_rst_val)[2] <- "length_mar"
dat_rst_val$type <- "Predicted"

### all fish recaptured in reach 1
dat_rc_val <- dat_rc[,c("tag_length", "rc_length", "tag_year")]
colnames(dat_rc_val)[2] <- "length_mar"
dat_rc_val$type <- "Observed"
dat_val <- rbind(dat_rst_val, dat_rc_val)


dat_17_val2 <- dat_rc[dat_rc$tag_year == 2017,]
dat_17_val2$pred_length <- dat_rst_17$pred_length[match(dat_17_val2$tag_length, dat_rst_17$tag_length)]
dat_17_val2$growth_var <- dat_rst_17$pred_growth_var[match(dat_17_val2$tag_length, dat_rst_17$tag_length)]

dat_18_val2 <- dat_rc[dat_rc$tag_year == 2018,]
dat_18_val2$pred_length <- dat_rst_18$pred_length[match(dat_18_val2$tag_length, dat_rst_18$tag_length)]
dat_18_val2$growth_var <- dat_rst_18$pred_growth_var[match(dat_18_val2$tag_length, dat_rst_18$tag_length)]

dat_19_val2 <- dat_rc[dat_rc$tag_year == 2019,]
dat_19_val2$pred_length <- dat_rst_19$pred_length[match(dat_19_val2$tag_length, dat_rst_19$tag_length)]
dat_19_val2$growth_var <- dat_rst_19$pred_growth_var[match(dat_19_val2$tag_length, dat_rst_19$tag_length)]

dat_val2 <- rbind(dat_17_val2, dat_18_val2, dat_19_val2)

# plot observed vs. predicted length at marine entry
ggplot(dat_val2) +
  coord_cartesian(xlim = c(65, 140), ylim = c(65, 140)) +
  geom_point(aes(x = rc_length, y = pred_length, color = tag_year)) +
  geom_segment(aes(x = rc_length, y = pred_length-growth_var, xend = rc_length, yend = pred_length+growth_var, color = tag_year)) +
  geom_abline(intercept = 0, slope = 1) +
  labs(x = "Observed length at marine entry (mm)", y = "Predicted length at marine entry (mm)", color = "Smolt year") +
  theme_classic()

# save validation plot
ggsave("fig_flmar_val.pdf", dpi = 300, height = 6, width = 6, units = "in")


