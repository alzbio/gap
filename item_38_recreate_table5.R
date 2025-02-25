###########################################################################
# Program Name: item_38_recreate_table5.R
# Description:  Recreate Table 5 from Mohs et al publication to find population used
# Client:       Global Alzheimers Platform (GAP)
# Protocol:     Bio-Hermes
#
# Author:       tduke
# Date Started: 2025-02-14
# Notes:
#
# Date Modified:
# Modified By:
# Reason Modified:
#
# Copyright 2025 Pentara Corporation
###########################################################################

# Attach packages ---------------------------------------------------------
library(tidyverse)
library(haven)
library(pROC)
library(pentarabeta)


# Read in the data --------------------------------------------------------
adlb = read_sas('P:/Clients/Global Alzheimers Platform (GAP)/Bio-Hermes/Data/Production/ADaM/adlb.sas7bdat')


# Define biomarkers and populations
biomarkers = c('TAU181P', 'TAU217P', 'BAB4240')
populations = c('Total', 'WHITEFL', 'HISPFL', 'BLACKFL')


results = data.frame()
for(i in biomarkers) {
  for(j in populations) {
    
    if(i == 'BAB4240') {
      data = adlb |>
        filter(LBNAM %ni% c('MERCK', 'QUANTERIX CORPORATION'))
    } else if(i == 'TAU181P') {
      data = adlb |>
        filter(LBNAM != 'ROCHE') 
    } else {
      data = adlb
    }
    
    # Filter and transform data
    data = data |>
      filter(ANL01FL == 'Y', AMYGR2 != '', 
             PARAMCD == i, DTYPE == '') |> # ANLSETFL == 'Y'
      select(USUBJID, AMYGR2, BLACKFL, WHITEFL, HISPFL, PARAMCD, AVAL) |>
      pivot_wider(names_from = PARAMCD, values_from = AVAL) |>
      mutate(AMYGR2 = if_else(AMYGR2 == 'Positive', 1, 0))
    
    ii = ensym(i)
    jj = ensym(j)
    
    # Filter by j if it isn't the Total Population, no filtering if it is
    if(j != 'Total') {
      data_to_use = data |>
        select(USUBJID, AMYGR2, !!ii, !!jj) |>
        filter(!!jj == 'Y') |>
        na.omit()
    } else {
      data_to_use = data |>
        select(USUBJID, AMYGR2, !!ii) |>
        na.omit()
    }
    
    n = length(unique(data_to_use$USUBJID))
    
    # Construct the model formula dynamically
    model_formula = reformulate(termlabels = i, response = 'AMYGR2')
    
    # Fit model and calculate AUC based on preds
    model = glm(model_formula, data = data, family = binomial)
    preds = predict(model, data_to_use, type = 'response')
    roc_obj = roc(response = data_to_use$AMYGR2, 
                  predictor = preds,
                  quiet = TRUE)
    
    # Extract values
    auc = round(auc(roc_obj), 4)
    ci = round(as.numeric(ci(roc_obj))[c(1, 3)], 4)
    se = round((ci[2] - ci[1]) / (2*qnorm(.975)), 4)
    
    result = data.frame(POP = j, 
                        BIO = i,
                        N = n,
                        AUC = auc, 
                        SE = se)
    results = rbind(results, result)
    
    # print(paste0('AUC for ', i, ' and ', j, ': ', round(auc(roc_obj), 4)))
    
  }
}

results |>
  mutate(POP = factor(POP, levels = c('Total', 'WHITEFL', 'HISPFL', 'BLACKFL'))) |>
  arrange(POP) |>
  mutate(IMP = 'No imputation',
         ADJUST = 'Unadjusted')
# These results match the non-imputed unadjusted values in Table 5


# Defining ANL01FL
anlfl = adlb |> 
  filter(PARAMCD %in% biomarkers, DTYPE == '', AMYGR2 != '') |> 
  select(USUBJID, PARAMCD, LBNAM, AMYGR2, ANL01FL, AVAL) |> 
  filter(case_when(
    PARAMCD == 'BAB4240' ~ LBNAM == 'C2N',
    PARAMCD == 'TAU181P' ~ LBNAM == 'QUANTERIX CORPORATION',
    PARAMCD == 'TAU217P' ~ LBNAM == 'LILLY CLINICAL DIAGNOSTICS LABORATORY'
  )) 

table(anlfl$ANL01FL)
