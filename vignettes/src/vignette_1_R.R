# Part I: Data Generation 
library(dplyr)
set.seed(123)
TT <- 20; GG <- 5
df <- data.frame(id = 1:(GG*TT))
df$G <- ((df$id-1) %% GG)+1
df$T <- floor((df$id-1)/GG)
df$id <- NULL
df <- df[order(df$G, df$T), ]

## One never treated group and 4 groups treated starting from period 4
## Treatment variable D:
#### Group 1 -> always 0
#### Group n -> 1 at period n+2

df$D <- 0
for (v in c(2,3,4,5)) {
    df$D <- ifelse(df$G == v & df$T == v+2, 1, df$D)
}
df <- df %>% group_by(.data$G) %>% mutate(D_stag = cumsum(.data$D)) %>% ungroup()

## The outcome is not-missing only every 4 years
df$Y <- ifelse(df$T %% 4 == 0, runif(n = nrow(df)) * (1 + 100*df$D_stag), NA)
df$D_stag <- NULL

## Let's browse the raw dataset
View(df)

## Now we have dataset with the issue of outcome being observed less frequently than the treatent

# Part II: Data Adjustment
## We need now to create the groups for the estimation.
## Keep in mind that, since the outcome is observable every 4 periods, we need to generate four indicators which point to (g,t) cells such that:
#### 1. treatment has never changed since the start of the panel or treatment has changed (at least once) at t and the change occurred in a non-missing outcome year
#### 2. treatment has never changed since the start of the panel or treatment has changed (at least once) at t and the change occurred the year before a non-missing outcome year
#### 3. treatment has never changed since the start of the panel or treatment has changed (at least once) at t and the change occurred two years before a non-missing outcome year
#### 4. treatment has never changed since the start of the panel or treatment has changed (at least once) at t and the change occurred three years before a non-missing outcome year

# a) Identify whether treatment has changed within each group (my suggestion: load the dplyr library)

df$D0 <- df$D[(df$G-1)*length(levels(factor(df$T)))+1] # We take the treatment value corresponding to T == 1 (If your groups are strings, just generate their id to make this method work)
df$D_change <- as.numeric(abs(df$D - df$D0) != 0)
df <- df %>% group_by(.data$G) %>% mutate(at_least_one_D_change = cumsum(.data$D_change)) %>% ungroup()
## As expected, the group 2 has at least one treatment change starting from period 4 (the period where its first treatment change occurred), group 3 from period 5 and so on.

# b) Identify when treatment changed for the first time

## We generate a variable (F_g) equal to the earliest period when there has been a treatment change
## For never switchers we use the same convention as de Chaisemartin & D'Haultfoeuille (2024) and set F_g = number of periods + 1
df <- df %>% group_by(.data$G) %>% 
mutate(never_treated = as.numeric(sum(.data$D_change, na.rm = TRUE) == 0)) %>%
mutate(F_g = ifelse(.data$never_treated == 1, max(df$T, na.rm = TRUE) +1, min(ifelse(.data$D_change == 0, NA, .data$T * .data$D_change), na.rm = TRUE))) %>% ungroup()

## Let's use the modulus operator to check whether the F_g falls on a 4th period or 1/2/3 years after (There are several ways to do this!)
## My approach (below): each year where the outcome is not missing takes value 1. Since the outcome is non-missing every 4 years, we just check whether the period is divisible by 4. If yes, the modulus operator by 4 yields 0. If not, it yield the remainder from dividing by 4. This remainder is always between 1 and 3. So, it is enough to increase the remainder by 1 to obtain 4 distinct values such that obs with non missing outcome (every fourth period) have value 1, obs after them value 2 and so on.
df$subsample <- (df$F_g %% 4) + 1

## From here, it is enough to multiply the subset variable by the indicator for at least one cange in D
df$subsample <- df$subsample * df$at_least_one_D_change


# Take a moment to look at the data and see yourself that the variable model_subset takes on values 0,1,2, 3 and 4. 
View(df)

# Keep only the observations with non-missing outcome
df <- subset(df, !is.na(df$Y))

## Part III: Estimation
## Notice that the estimation will be performed with the staggerized version of the treatment (that we called at_least_one_D_change)

library(DIDmultiplegtDYN)
effects <- 2
table <- NULL
for (j in 1:4) {
    temp <- did_multiplegt_dyn(
        subset(df, df$subsample %in% c(0, j)), "Y", "G", "T", "at_least_one_D_change", 
        graph_off = TRUE, effects = effects)
    rownames(temp$results$Effects) <- 
        sapply(1:temp$results$N_Effects, function(x) paste0("Effect_",  j + (x-1) * 4))
    table <- rbind(table, temp$results$Effects)
}
rown <- unlist(strsplit(rownames(table), "_")) 
table <- cbind(table, as.numeric(rown[rown != "Effect"]))
print(table[order(table[,ncol(table)]),1:(ncol(table)-1)])

# Bigger Data for graphs
set.seed(123)
library(dplyr)
library(DIDmultiplegtDYN)

TT <- 20; GG <- 1000
df <- data.frame(id = 1:(GG*TT))
df$G <- ((df$id-1) %% GG)+1
df$T <- floor((df$id-1)/GG)
df$id <- NULL
df <- df[order(df$G, df$T), ]
df$D <- 0
df$D <- ifelse(df$T == ((df$G - 1) %% 5) + 3  & (df$G - 1) %% 5 != 0, 1, df$D)
df <- df %>% group_by(.data$G) %>% mutate(D_stag_temp = cumsum(.data$D)) %>% ungroup()
df <- df %>% group_by(.data$G) %>% mutate(D_stag = cumsum(.data$D_stag_temp)) %>% ungroup()
df$Y <- ifelse(df$T %% 4 == 0, runif(n = nrow(df)) * (1 + df$D_stag), NA)
df$D_stag <- NULL
df$D0 <- df$D[(df$G-1)*length(levels(factor(df$T)))+1]
df$D_change <- as.numeric(abs(df$D - df$D0) != 0)
df <- df %>% group_by(.data$G) %>% mutate(at_least_one_D_change = cumsum(.data$D_change)) %>% ungroup()
df <- df %>% group_by(.data$G) %>% 
mutate(never_treated = as.numeric(sum(.data$D_change, na.rm = TRUE) == 0)) %>%
mutate(F_g = ifelse(.data$never_treated == 1, max(df$T, na.rm = TRUE) +1, min(ifelse(.data$D_change == 0, NA, .data$T * .data$D_change), na.rm = TRUE))) %>% ungroup()
df$subsample <- (4 - (df$F_g %% 4)) * (df$F_g %% 4 != 0) + 1
df$subsample <- df$subsample * df$at_least_one_D_change
df <- subset(df, !is.na(df$Y))

effects <- 2
table <- NULL
for (j in 1:4) {
    temp <- did_multiplegt_dyn(
        subset(df, df$subsample %in% c(0, j)), "Y", "G", "T", "at_least_one_D_change", 
        graph_off = TRUE, effects = effects)
    rownames(temp$results$Effects) <- 
        sapply(1:temp$results$N_Effects, function(x) paste0("Effect_",  j + (x-1) * 4))
    table <- rbind(table, temp$results$Effects)
}
rown <- unlist(strsplit(rownames(table), "_")) 
table <- cbind(table, as.numeric(rown[rown != "Effect"]))
print(table[order(table[,ncol(table)]),1:(ncol(table)-1)])

library(ggplot2)
table <- table[order(table[,ncol(table)]), ]
table <- rbind(rep(0, ncol(table)), table)
colnames(table)[ncol(table)] <- "Time"
table <- as.data.frame(table)
out_plot <- ggplot(table, aes(x = .data$Time, y = .data$Estimate, group = 1)) + 
geom_line(colour = "blue") +
geom_errorbar(data = ~dplyr::filter(.x, table$Estimate != 0), aes(ymin = .data[["LB CI"]], ymax = .data[["UB CI"]]), 
position=position_dodge(0.05), width = 0.2, colour = "red") + 
geom_point(colour = "blue") + 
ggtitle("DID, from last period before treatment changes (t=0) to t") + 
xlab("Relative time to last period before treatment changes (t=0)") +
theme(plot.title = element_text(hjust = 0.5))
print(out_plot)

# Bigger Data for graphs - PLACEBO
set.seed(123)
library(dplyr)
library(DIDmultiplegtDYN)

TT <- 20; GG <- 1000
df <- data.frame(id = 1:(GG*TT))
df$G <- ((df$id-1) %% GG)+1
df$T <- floor((df$id-1)/GG)
df$id <- NULL
df <- df[order(df$G, df$T), ]
df$D <- 0
df$D <- ifelse(df$T == ((df$G - 1) %% 5) + 3  & (df$G - 1) %% 5 != 0, 1, df$D)
df <- df %>% group_by(.data$G) %>% mutate(D_stag_temp = cumsum(.data$D)) %>% ungroup()
df <- df %>% group_by(.data$G) %>% mutate(D_stag = cumsum(.data$D_stag_temp)) %>% ungroup()
df$Y <- ifelse(df$T %% 4 == 0, runif(n = nrow(df)) * (1 + df$D_stag), NA)
df$D_stag <- NULL
df$D0 <- df$D[(df$G-1)*length(levels(factor(df$T)))+1]
df$D_change <- as.numeric(abs(df$D - df$D0) != 0)
df <- df %>% group_by(.data$G) %>% mutate(at_least_one_D_change = cumsum(.data$D_change)) %>% ungroup()
df <- df %>% group_by(.data$G) %>% 
mutate(never_treated = as.numeric(sum(.data$D_change, na.rm = TRUE) == 0)) %>%
mutate(F_g = ifelse(.data$never_treated == 1, max(df$T, na.rm = TRUE) +1, min(ifelse(.data$D_change == 0, NA, .data$T * .data$D_change), na.rm = TRUE))) %>% ungroup()
df$subsample <- (4 - (df$F_g %% 4)) * (df$F_g %% 4 != 0) + 1
df$subsample <- df$subsample * df$at_least_one_D_change
df <- subset(df, !is.na(df$Y))

effects <- 2
table <- NULL
for (j in 1:4) {
    temp <- did_multiplegt_dyn(
        subset(df, df$subsample %in% c(0, j)), "Y", "G", "T", "at_least_one_D_change", 
        graph_off = TRUE, effects = effects)
    rownames(temp$results$Effects) <- 
        sapply(1:temp$results$N_Effects, function(x) paste0("Effect_",  j + (x-1) * 4))
    table <- rbind(table, temp$results$Effects)
}
rown <- unlist(strsplit(rownames(table), "_")) 
table <- cbind(table, as.numeric(rown[rown != "Effect"]))
print(table[order(table[,ncol(table)]),1:(ncol(table)-1)])

library(ggplot2)
table <- table[order(table[,ncol(table)]), ]
table <- rbind(rep(0, ncol(table)), table)
colnames(table)[ncol(table)] <- "Time"
table <- as.data.frame(table)
out_plot <- ggplot(table, aes(x = .data$Time, y = .data$Estimate, group = 1)) + 
geom_line(colour = "blue") +
geom_errorbar(data = ~dplyr::filter(.x, table$Estimate != 0), aes(ymin = .data[["LB CI"]], ymax = .data[["UB CI"]]), 
position=position_dodge(0.05), width = 0.2, colour = "red") + 
geom_point(colour = "blue") + 
ggtitle("DID, from last period before treatment changes (t=0) to t") + 
xlab("Relative time to last period before treatment changes (t=0)") +
theme(plot.title = element_text(hjust = 0.5))
print(out_plot)