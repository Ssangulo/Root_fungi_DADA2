### NEW beta regression code integrating all orders

library(betareg)
library(dplyr)
library(tidyr)
library(ggplot2)

ps <- individual_ps

### ============================================================
### 1. Extract OTU table (rows = OTUs, columns = samples)
### ============================================================

otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) {
  otu <- t(otu)
}

tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)

# Safe extraction of Order column
if ("Order" %in% names(tax)) {
  order_vec <- tolower(tax$Order)
  order_vec[is.na(order_vec)] <- ""
} else {
  stop("No 'Order' column found in tax table.")
}

orders_target <- c(
  "helotiales",
  "sebacinales",
  "chaetothyriales",
  "sclerococcales",
  "agaricales",
  "pleosporales",
  "hypocreales"
)

get_otus_for_order <- function(ord) {
  rownames(tax)[grepl(ord, order_vec)]
}

order_otu_list <- lapply(orders_target, get_otus_for_order)
names(order_otu_list) <- orders_target

# sanity check
sapply(order_otu_list, length)

### ============================================================
### 3. Sample metadata aligned to OTU table
### ============================================================

sdat    <- as.data.frame(sample_data(ps))
samples <- colnames(otu)

total_counts <- colSums(otu)


### ============================================================
### 4. Helper to build proper beta-regression dataframes
### ============================================================

make_beta_df <- function(order_name, otu_table, otu_ids, total_counts, sdat, samples) {
  
  if (length(otu_ids) == 0) {
    warning(paste("Order has zero OTUs:", order_name))
    return(NULL)
  }
  
  counts <- colSums(otu_table[otu_ids, , drop = FALSE])
  
  df <- data.frame(
    sample    = samples,
    count     = counts[samples],
    total     = total_counts[samples],
    site      = sdat[samples, "site", drop = TRUE],
    elevation = sdat[samples, "elevation", drop = TRUE],
    stringsAsFactors = FALSE
  )
  
  # raw relative abundance
  df$prop_raw <- ifelse(df$total > 0, df$count / df$total, NA_real_)
  
  # Smithson-Verkuilen adjustment (0,1)
  n <- nrow(df)
  df$Proportion <- (df$prop_raw * (n - 1) + 0.5) / n
  
  df$site      <- factor(df$site)
  df$elevation <- as.numeric(df$elevation)
  
  df
}


### ============================================================
### 5. Build ALL beta-regression dataframes
### ============================================================

df_list <- list()

for (ord in orders_target) {
  otu_ids <- order_otu_list[[ord]]
  df_list[[ord]] <- make_beta_df(ord, otu, otu_ids, total_counts, sdat, samples)
}

### ============================================================
### 6. Fit beta regression models for each order
### ============================================================

model_list <- list()

for (ord in orders_target) {
  df_ord <- df_list[[ord]]
  if (!is.null(df_ord)) {
    model_list[[ord]] <- betareg(Proportion ~ elevation + site, data = df_ord)
  } else {
    model_list[[ord]] <- NULL
  }
}

# individual models
mod_helot  <- model_list$helotiales
mod_seba   <- model_list$sebacinales
mod_chaet  <- model_list$chaetosphaeriales
mod_sclero <- model_list$sclerococcales
mod_agar   <- model_list$agaricales
mod_pleo   <- model_list$pleosporales
mod_hypo   <- model_list$hypocreales

# usual usage:
summary(mod_seba)
coef(mod_seba)

# function to extract elevation coefficient from a betareg model
extract_elev_coef <- function(mod, order_name) {
  if (is.null(mod)) return(NULL)
  
  sm <- summary(mod)
  coefs <- sm$coefficients$mean  # mean submodel (not precision)
  
  if (!("elevation" %in% rownames(coefs))) {
    warning(paste("No 'elevation' term in model for order:", order_name))
    return(NULL)
  }
  
  elev <- coefs["elevation", ]
  
  data.frame(
    Order      = order_name,
    term       = "elevation",
    estimate   = elev["Estimate"],
    se         = elev["Std. Error"],
    z          = elev["z value"],
    p_value    = elev["Pr(>|z|)"],
    pseudo_R2  = sm$pseudo.r.squared,
    stringsAsFactors = FALSE
  )
}

# loop over all models and bind into one table
elev_summary <- bind_rows(
  lapply(names(model_list), function(ord) {
    extract_elev_coef(model_list[[ord]], ord)
  })
)

elev_summary

write_xlsx(
  elev_summary,
  "/data/lastexpansion/danieang/models/elevation_beta_orders.xlsx"
)

sm <- summary(mod_seba)
mean_coefs <- as.data.frame(sm$coefficients$mean)

write_xlsx(
  mean_coefs,
  "/data/lastexpansion/danieang/models/betareg_mod_seba.xlsx"
)
