--- 
author: "Ushani Atapattu"
date: "`r Sys.Date()`"
editor_options: 
  markdown: 
    wrap: 80
output: 
  code_folding: show
  html_document: 
  number_sections: true
  toc: true
  toc_float: true
title: "Molecular diagnostic survey of the prevalence and associated risk factors for tick-borne pathogens in pet dogs in Sri Lanka"

---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r include=FALSE}
#Load data
ticks_data <- read.csv("tbd_data.csv")

```

## Libraries

```{r include=FALSE}
library(tidyverse)
library(lubridate)
library(janitor)
library(skimr)
#library(epiR)
library(epitools)
#library(binom)
library(ggpubr)
library(MASS)
library(sjPlot)
library(flextable)
#library(UpSetR)
library(broom)
library(lme4)
```

## Data cleaning

```{r echo=FALSE}

ticks_clean <- ticks_data %>% 
  mutate(across(c("dog_age_mo", "dog_breed"), ~if_else(.=="unknwn", NA_character_, .))) %>%
  mutate(across(c( #Missing values
      "dog_sex",
      "dog_neuterstat",
      "dog_breedgrp",
      "dog_agegrp",
      "tk_roa",
      "tk_systemic",
      "tk_topical",
      "sk_tick",
      "sk_flea",
      "sk_lice",
      "dog_bcs",
      "dog_condition",
      "dog_tkrx",
      "dog_dewrmed",
      "cs_pyrexia",
      "cs_mucsmem",
      "cs_ln",
      "cs_breath",
      "cs_sphep",
      "cs_appetite",
      "cs_ascites",
      "sk_dx",
      "hx_tkfvr",
      "cs_ears"
    ), ~if_else(.=="na", NA_character_, .))) %>%
    mutate( dog_breedgrp2 = case_when( # Breed groups as local and exotic
    dog_breedgrp == 1 | dog_breedgrp == 3 ~ 1,
    dog_breedgrp == 2 | dog_breedgrp == 4 ~ 2,
     TRUE ~ NA_real_  ))%>%
  mutate(dog_breedgrp = replace(dog_breedgrp, dog_breedgrp %in% c(3, 4), 1))%>%
  mutate(dog_age_mo = as.numeric(dog_age_mo)) %>%
  mutate(tbp_infection = ifelse(tbp_infecnum >= 1, 1, 0)) %>%
  mutate(across(c( # converting into "factor" data class
    "tk_roa", 
    "hx_tkfvr",
    "dog_tkrx", 
    "cs_pyrexia",
    "dog_province",
    "tbp_anaplasma",
    "tbp_ecanis",
    "tbp_mycoplasma",
    "dog_breedgrp2",
    "tbp_bgibsoni",
    "tbp_bvogeli",
    "tbp_hcanis",
    "dog_zone",
    "tbp_infection",
    "dog_sex",
    "dog_neuterstat",
    "dog_bcs",
    "dog_condition",
    "dog_zone2",
    "dog_tkrx",
    "tbp_infecnum",
    "tbp_infstat"
  ), as_factor)) %>%
  mutate(
    dog_sexstat = case_when(
      dog_sex == 1 & dog_neuterstat == 1 ~ 1,
      dog_sex == 1 & dog_neuterstat == 0 ~ 2,
      dog_sex == 0 & dog_neuterstat == 1 ~ 3,
      dog_sex == 0 & dog_neuterstat == 0 ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(dog_sexstat = as_factor(dog_sexstat)) %>%
  mutate(
    cs_tickfever = case_when(
      cs_pyrexia == 1 | cs_mucsmem == 1 | cs_ln == 1 | cs_sphep == 1 | cs_appetite == 1 ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(cs_tickfever = as_factor(cs_tickfever))%>%
  mutate(bdt_pos = case_when( 
    tbp_anaplasma == 1 | tbp_bvogeli == 1 | tbp_hcanis == 1 | tbp_ecanis == 1 ~ 1,
    TRUE ~ 0)
    ) %>%
  mutate(bdt_pos = as_factor(bdt_pos)) %>%
  mutate(tbp_fight = case_when(tbp_bgibsoni == 1 | tbp_mycoplasma == 1 ~ 1,
                               TRUE ~ 0)) %>%
   mutate(tbp_fight = as_factor(tbp_fight)) %>%
  mutate(dog_ageadult = case_when(dog_agegrp == 3 | dog_agegrp == 4 | dog_agegrp == 5 | dog_agegrp == 6 ~ 1,
                               TRUE ~ NA_real_)) %>%
   mutate(dog_ageadult = as_factor(dog_ageadult))%>%
  mutate(cs_abmucsmemb = case_when( cs_mucsmem == 1 | cs_mucsmem == 2 ~ 1,
                                     cs_mucsmem  == 3~ 0,
                                     TRUE ~ NA_real_))%>%
  mutate(cs_abtemp = case_when( cs_pyrexia == 1 | cs_pyrexia == 2 ~ 1,
                                cs_pyrexia == 3 ~ 0,
                                     TRUE ~ NA_real_))%>%
    mutate(cs_lowbcs = case_when(dog_bcs == 1 | dog_bcs == 2 ~ 1,
                                 dog_bcs == 3| dog_bcs == 4 ~ 0,
                                     TRUE ~ NA_real_))%>%
  mutate(tbp_mixed = case_when(tbp_infecnum == 1 ~ 1,
                               tbp_infecnum == 2 | tbp_infecnum == 3 |tbp_infecnum == 4 ~ 2,
                               tbp_infecnum == 0 ~ 0,
                               TRUE ~ NA_real_))%>%
  mutate(cs_mixed = as_factor(tbp_mixed))%>%
   mutate(tbp_infecnumgrp = case_when(
    tbp_infecnum == "0" ~ "0", #no infection
    tbp_infecnum == "1" ~ "1",# mono-infections
    tbp_infecnum == "2" ~ "2",# mixed infections
    tbp_infecnum == "3" ~ "2",
    tbp_infecnum == "4" ~ "2"
     ))%>%
  mutate(dog_age_yr = dog_age_mo/12)%>%
  mutate(tbp_infecnumgrp = as_factor(tbp_infecnumgrp))
# geo-climatic zones : 1=mid/up country wet zone; 2=low country wet zone; 3=low country dry zone
```

```{r}

ticks_clean$tbp_infecnumgrp <- factor(ticks_clean$tbp_infecnumgrp, levels = c("0","1","2"))
ticks_clean$tbp_infecnumgrp <- ordered(ticks_clean$tbp_infecnumgrp, levels = c("0","1","2"))

```

```{r}

#Function to tabulate data in a table

tab_data <- function(df, x, y) {
  tab <-table(df[,x], df[,y], useNA = "ifany")
  print(tab)
}
```

## Sample descriptives according to geo-climatic zone

```{r}


desc_age <- tab_data(ticks_clean, y ="dog_zone2", x = "dog_agegrp")
desc_sex <- tab_data(ticks_clean, y ="dog_zone2", x = "dog_sex")
desc_breed <- tab_data(ticks_clean, y ="dog_zone2", x = "dog_breedgrp2")
desc_neuterstat <- tab_data(ticks_clean, y ="dog_zone2", x = "dog_neuterstat")
desc_tick <- tab_data(ticks_clean, y ="dog_zone2", x = "sk_tick")
desc_flea <- tab_data(ticks_clean, y ="dog_zone2", x = "sk_flea")
desc_tickrx <- tab_data(ticks_clean, y ="dog_zone2", x = "dog_tkrx")

```

## TBP infection in the study group

```{r}
# Function to tabulate and calculate the confidence intervals

sum_data <- function(df, x, y){
  
  # tabulate data in a contigency table
  tab <-table(df[ ,x], df[ ,y], useNA = "ifany")
  
  #calculate the total from each group
  
 total <- marginSums(tab, margin = 2)
 
 
 tab <- rbind(tab, total)

#calculate propotion infected with 95% confidence intervals 
 
output <- mapply(binom.exact, tab[2, ], tab[3, ], conf.level = 0.95)

print(output)

}

```

### Overall TBP infection

```{r}

prop_age <-sum_data(ticks_clean, x ="tbp_infstat", y = "dog_agegrp")
prop_sex <- sum_data(ticks_clean, x ="tbp_infstat", y = "dog_sex")
prop_breed <- sum_data(ticks_clean, x ="tbp_infstat", y = "dog_breedgrp2")
prop_neuterstat <- sum_data(ticks_clean, x ="tbp_infstat", y = "dog_neuterstat")
prop_tick <- sum_data(ticks_clean, x ="tbp_infstat", y = "sk_tick")
prop_flea <- sum_data(ticks_clean, x="tbp_infstat", y = "sk_flea")
prop_tickrx <- sum_data(ticks_clean, x ="tbp_infstat", y = "dog_tkrx")

```

### *A. platys* infection

```{r}
prop_age <-sum_data(ticks_clean, x ="tbp_anaplasma", y = "dog_agegrp")
prop_sex <- sum_data(ticks_clean, x ="tbp_anaplasma", y = "dog_sex")
prop_breed <- sum_data(ticks_clean, x ="tbp_anaplasma", y = "dog_breedgrp2")
prop_neuterstat <- sum_data(ticks_clean, x ="tbp_anaplasma", y = "dog_neuterstat")
prop_tick <- sum_data(ticks_clean, x ="tbp_anaplasma", y = "sk_tick")
prop_flea <- sum_data(ticks_clean, x="tbp_anaplasma", y = "sk_flea")
prop_tickrx <- sum_data(ticks_clean, x ="tbp_anaplasma", y = "dog_tkrx")
```

### *E. canis* infection

```{r}
prop_age <-sum_data(ticks_clean, x ="tbp_ecanis", y = "dog_agegrp")
prop_sex <- sum_data(ticks_clean, x ="tbp_ecanis", y = "dog_sex")
prop_breed <- sum_data(ticks_clean, x ="tbp_ecanis", y = "dog_breedgrp2")
prop_neuterstat <- sum_data(ticks_clean, x ="tbp_ecanis", y = "dog_neuterstat")
prop_tick <- sum_data(ticks_clean, x ="tbp_ecanis", y = "sk_tick")
prop_flea <- sum_data(ticks_clean, x="tbp_ecanis", y = "sk_flea")
prop_tickrx <- sum_data(ticks_clean, x ="tbp_ecanis", y = "dog_tkrx")
```

### Hemotropic mycoplasma infection

```{r}
prop_age <-sum_data(ticks_clean, x ="tbp_mycoplasma", y = "dog_agegrp")
prop_sex <- sum_data(ticks_clean, x ="tbp_mycoplasma", y = "dog_sex")
prop_breed <- sum_data(ticks_clean, x ="tbp_mycoplasma", y = "dog_breedgrp2")
prop_neuterstat <- sum_data(ticks_clean, x ="tbp_mycoplasma", y = "dog_neuterstat")
prop_tick <- sum_data(ticks_clean, x ="tbp_mycoplasma", y = "sk_tick")
prop_flea <- sum_data(ticks_clean, x="tbp_mycoplasma", y = "sk_flea")
prop_tickrx <- sum_data(ticks_clean, x ="tbp_mycoplasma", y = "dog_tkrx")
```

### *B. gibsoni* infection

```{r}
prop_age <-sum_data(ticks_clean, x ="tbp_bgibsoni", y = "dog_agegrp")
prop_sex <- sum_data(ticks_clean, x ="tbp_bgibsoni", y = "dog_sex")
prop_breed <- sum_data(ticks_clean, x ="tbp_bgibsoni", y = "dog_breedgrp2")
prop_neuterstat <- sum_data(ticks_clean, x ="tbp_bgibsoni", y = "dog_neuterstat")
prop_tick <- sum_data(ticks_clean, x ="tbp_bgibsoni", y = "sk_tick")
prop_flea <- sum_data(ticks_clean, x="tbp_bgibsoni", y = "sk_flea")
prop_tickrx <- sum_data(ticks_clean, x ="tbp_bgibsoni", y = "dog_tkrx")
```

### *B. vogeli* infection

```{r}
prop_age <-sum_data(ticks_clean, x ="tbp_bvogeli", y = "dog_agegrp")
prop_sex <- sum_data(ticks_clean, x ="tbp_bvogeli", y = "dog_sex")
prop_breed <- sum_data(ticks_clean, x ="tbp_bvogeli", y = "dog_breedgrp2")
prop_neuterstat <- sum_data(ticks_clean, x ="tbp_bvogeli", y = "dog_neuterstat")
prop_tick <- sum_data(ticks_clean, x ="tbp_bvogeli", y = "sk_tick")
prop_flea <- sum_data(ticks_clean, x="tbp_bvogeli", y = "sk_flea")
prop_tickrx <- sum_data(ticks_clean, x ="tbp_bvogeli", y = "dog_tkrx")
```

### *H. canis* infection

```{r}
prop_age <-sum_data(ticks_clean, x ="tbp_hcanis", y = "dog_agegrp")
prop_sex <- sum_data(ticks_clean, x ="tbp_hcanis", y = "dog_sex")
prop_breed <- sum_data(ticks_clean, x ="tbp_hcanis", y = "dog_breedgrp2")
prop_neuterstat <- sum_data(ticks_clean, x ="tbp_hcanis", y = "dog_neuterstat")
prop_tick <- sum_data(ticks_clean, x ="tbp_hcanis", y = "sk_tick")
prop_flea <- sum_data(ticks_clean, x="tbp_hcanis", y = "sk_flea")
prop_tickrx <- sum_data(ticks_clean, x ="tbp_hcanis", y = "dog_tkrx")
```

### Mixed infections

```{r}

#tabulate mixed infection propotions

sum_data_type <- function(df, x, y){
  
  # tabulate data in a contigency table
  tab <-table(df[ ,x], df[ ,y], useNA = "ifany")
  
  
  #callculate the total from each group
  
 total <- marginSums(tab, margin = 1)
  
 
 tab <- cbind(tab, total)
 
 
 #print(tab)

#calculate propotion infected with 95% confidence intervals 
 # Single infection
single <- mapply(binom.exact, tab[ ,1], tab[ ,4], conf.level = 0.95)
 # Mixed infection
mixed <- mapply(binom.exact, tab[ ,2], tab[ ,4], conf.level = 0.95)
print("single infection", quotes = FALSE)
print(single)
print("mixed infection", quotes = FALSE)
print(mixed)

}

```

```{r}

prop_age_tbp <- sum_data_type(df=ticks_clean, x="dog_agegrp", y="tbp_infecnumgrp")
prop_sex_tbp <- sum_data_type(df=ticks_clean, x="dog_sex", y="tbp_infecnumgrp")
prop_breed_tbp <- sum_data_type(df=ticks_clean, x="dog_breedgrp2", y="tbp_infecnumgrp")
prop_neuterstat_tbp <- sum_data_type(df=ticks_clean, x="dog_neuterstat", y="tbp_infecnumgrp")
prop_tick_tbp <- sum_data_type(df=ticks_clean, x="sk_tick", y="tbp_infecnumgrp")
prop_flea_tbp <- sum_data_type(df=ticks_clean, x="sk_flea", y="tbp_infecnumgrp")
prop_tkrx_tbp <- sum_data_type(df=ticks_clean, x="dog_tkrx", y="tbp_infecnumgrp")
prop_zone_tbp <- sum_data_type(df=ticks_clean, x="dog_zone2", y="tbp_infecnumgrp")
```

## Univariable associations

```{r}
# Code for univariable association with fishers exact test or chi square test with significance at 0.2

test_uni     <-  function(df, x, y, alpha = 0.2) {
  
  #tabulate data in a contigency table
  
  tab             <- table(df[, x], df[, y])
  
 # print(tab)
  
  total_cells     <- nrow(tab) * ncol(tab)
  
  cst             <- chisq.test(tab, correct = FALSE)
  
  expected_counts <- cst$expected

  cells_gte_five  <- sum(expected_counts >= 5)
  
  print("Odds Ratio", quotes = FALSE)
 odds             <- epitools::oddsratio(x=tab, method="wald") %>%
        print()
 or               <- odds$measure
variable     <- c(x,x)
  
  # if the assumptions of the chi squared test are violated...
  if (((cells_gte_five / total_cells)  < 0.8) || any(expected_counts < 1)) {
    
    fet           <- fisher.test(tab)
 

   p         <- fet[["p.value"]]
   
   selected <- ifelse(fet[["p.value"]]<=alpha, "yes", "no")
 
      or  <- cbind(variable, or, p, selected)%>%
      as.data.frame()%>%
     print()
    
      print("Fischer's test used", quotes = FALSE)
  } else {
    
     p  <- cst[["p.value"]]
     selected <- ifelse(cst[["p.value"]]<=alpha, "yes", "no")
   
      or <- cbind(variable, or, p, selected )%>%
      as.data.frame()%>%
      print()
      
      print("chi square test used", quotes = FALSE)
      }


}

```

1.  **Host variables:**

-   Age

-   Breed

-   Sex

-   Neuter status

-   Tick infestation

-   Flea infestation

2.  **Geo-climatic zone**
3.  **Use of ectoparasiticides**

### *A. platys*

```{r}


uni_age_ana <- test_uni(df=ticks_clean, x="dog_agegrp", y="tbp_anaplasma")
uni_sex_ana <- test_uni(df=ticks_clean, x="dog_sex", y="tbp_anaplasma")
uni_breed_ana <- test_uni(df=ticks_clean, x="dog_breedgrp2", y="tbp_anaplasma")
uni_neuterstat_ana <- test_uni(df=ticks_clean, x="dog_neuterstat", y="tbp_anaplasma")
uni_tick_ana <- test_uni(df=ticks_clean, x="sk_tick", y="tbp_anaplasma")
uni_flea_ana <- test_uni(df=ticks_clean, x="sk_flea", y="tbp_anaplasma")
uni_tkrx_ana <- test_uni(df=ticks_clean, x="dog_tkrx", y="tbp_anaplasma")
uni_zone_ana <- test_uni(df=ticks_clean, x="dog_zone2", y="tbp_anaplasma")


```

### *E. canis*

```{r}

uni_age_ec <- test_uni(df=ticks_clean, x="dog_agegrp", y="tbp_ecanis")
uni_sex_ec <- test_uni(df=ticks_clean, x="dog_sex", y="tbp_ecanis")
uni_breed_ec <- test_uni(df=ticks_clean, x="dog_breedgrp2", y="tbp_ecanis")
uni_neuterstat_ec <- test_uni(df=ticks_clean, x="dog_neuterstat", y="tbp_ecanis")
uni_tick_ec <- test_uni(df=ticks_clean, x="sk_tick", y="tbp_ecanis")
uni_flea_ec <- test_uni(df=ticks_clean, x="sk_flea", y="tbp_ecanis")
uni_tkrx_ec <- test_uni(df=ticks_clean, x="dog_tkrx", y="tbp_ecanis")
uni_zone_ec <- test_uni(df=ticks_clean, x="dog_zone2", y="tbp_ecanis")
```

### Hemotropic mycoplasma spp.

```{r}

uni_age_myc <- test_uni(df=ticks_clean, x="dog_agegrp", y="tbp_mycoplasma")
uni_sex_myc <- test_uni(df=ticks_clean, x="dog_sex", y="tbp_mycoplasma")
uni_breed_myc <- test_uni(df=ticks_clean, x="dog_breedgrp2", y="tbp_mycoplasma")
uni_neuterstat_myc <- test_uni(df=ticks_clean, x="dog_neuterstat", y="tbp_mycoplasma")
uni_tick_myc <- test_uni(df=ticks_clean, x="sk_tick", y="tbp_mycoplasma")
uni_flea_myc <- test_uni(df=ticks_clean, x="sk_flea", y="tbp_mycoplasma")
uni_tkrx_myc <- test_uni(df=ticks_clean, x="dog_tkrx", y="tbp_mycoplasma")
uni_zone_myc <- test_uni(df=ticks_clean, x="dog_zone2", y="tbp_mycoplasma")

```

### *B. gibsoni*

```{r}

uni_age_bg <- test_uni(df=ticks_clean, x="dog_agegrp", y="tbp_bgibsoni")
uni_sex_bg <- test_uni(df=ticks_clean, x="dog_sex", y="tbp_bgibsoni")
uni_breed_bg <- test_uni(df=ticks_clean, x="dog_breedgrp2", y="tbp_bgibsoni")
uni_neuterstat_bg <- test_uni(df=ticks_clean, x="dog_neuterstat", y="tbp_bgibsoni")
uni_tick_bg <- test_uni(df=ticks_clean, x="sk_tick", y="tbp_bgibsoni")
uni_flea_bg <- test_uni(df=ticks_clean, x="sk_flea", y="tbp_bgibsoni")
uni_tkrx_bg <- test_uni(df=ticks_clean, x="dog_tkrx", y="tbp_bgibsoni")
uni_zone_bg <- test_uni(df=ticks_clean, x="dog_zone2", y="tbp_bgibsoni")
```

### *B. vogeli*

```{r}


uni_age_bv <- test_uni(df=ticks_clean, x="dog_agegrp", y="tbp_bvogeli")
uni_sex_bv <- test_uni(df=ticks_clean, x="dog_sex", y="tbp_bvogeli")
uni_breed_bv <- test_uni(df=ticks_clean, x="dog_breedgrp2", y="tbp_bvogeli")
uni_neuterstat_bv <- test_uni(df=ticks_clean, x="dog_neuterstat", y="tbp_bvogeli")
uni_tick_bv <- test_uni(df=ticks_clean, x="sk_tick", y="tbp_bvogeli")
uni_flea_bv <- test_uni(df=ticks_clean, x="sk_flea", y="tbp_bvogeli")
uni_tkrx_bv <- test_uni(df=ticks_clean, x="dog_tkrx", y="tbp_bvogeli")
uni_zone_bv <- test_uni(df=ticks_clean, x="dog_zone2", y="tbp_bvogeli")


```

### *H. canis*

```{r}

uni_age_hc <- test_uni(df=ticks_clean, x="dog_agegrp", y="tbp_hcanis")
uni_sex_hc <- test_uni(df=ticks_clean, x="dog_sex", y="tbp_hcanis")
uni_breed_hc <- test_uni(df=ticks_clean, x="dog_breedgrp2", y="tbp_hcanis")
uni_neuterstat_hc <- test_uni(df=ticks_clean, x="dog_neuterstat", y="tbp_hcanis")
uni_tick_hc <- test_uni(df=ticks_clean, x="sk_tick", y="tbp_hcanis")
uni_flea_hc <- test_uni(df=ticks_clean, x="sk_flea", y="tbp_hcanis")
uni_tkrx_hc <- test_uni(df=ticks_clean, x="dog_tkrx", y="tbp_hcanis")
uni_zone_hc <- test_uni(df=ticks_clean, x="dog_zone2", y="tbp_hcanis")

```

### Overall TBP infection

```{r}

uni_age_tbp <- test_uni(df=ticks_clean, x="dog_agegrp", y="tbp_infstat")
uni_sex_tbp <- test_uni(df=ticks_clean, x="dog_sex", y="tbp_infstat")
uni_breed_tbp <- test_uni(df=ticks_clean, x="dog_breedgrp2", y="tbp_infstat")
uni_neuterstat_tbp <- test_uni(df=ticks_clean, x="dog_neuterstat", y="tbp_infstat")
uni_tick_tbp <- test_uni(df=ticks_clean, x="sk_tick", y="tbp_infstat")
uni_flea_tbp <- test_uni(df=ticks_clean, x="sk_flea", y="tbp_infstat")
uni_tkrx_tbp <- test_uni(df=ticks_clean, x="dog_tkrx", y="tbp_infstat")
uni_zone_tbp <- test_uni(df=ticks_clean, x="dog_zone2", y="tbp_infstat")

```

### Mixed infection

```{r}

uni_age_infnum <- test_uni(df=ticks_clean, x="dog_agegrp", y="tbp_infecnumgrp")
uni_sex_infnum <- test_uni(df=ticks_clean, x="dog_sex", y="tbp_infecnumgrp")
uni_breed_infnum <- test_uni(df=ticks_clean, x="dog_breedgrp2", y="tbp_infecnumgrp")
uni_neuterstat_infnum <- test_uni(df=ticks_clean, x="dog_neuterstat", y="tbp_infecnumgrp")
uni_tick_infnum <- test_uni(df=ticks_clean, x="sk_tick", y="tbp_infecnumgrp")
uni_flea_infnum <- test_uni(df=ticks_clean, x="sk_flea", y="tbp_infecnumgrp")
uni_tkrx_infnum <- test_uni(df=ticks_clean, x="dog_tkrx", y="tbp_infecnumgrp")
uni_zone_infnum <- test_uni(df=ticks_clean, x="dog_zone2", y="tbp_infecnumgrp")

```

## Clinical Signs and TBP infection - Univariable Analysis

```{r}
# Code for univariable association with fishers exact test or chi square test with significance at 0.2

test_uni_cs     <-  function(df, x, y, alpha = 0.05) {
  
  #tabulate data in a contigency table
  
  tab             <- table(df[, x], df[, y])
  
 # print(tab)
  
  total_cells     <- nrow(tab) * ncol(tab)
  
  cst             <- chisq.test(tab, correct = FALSE)
  
  expected_counts <- cst$expected

  cells_gte_five  <- sum(expected_counts >= 5)
  
  print("Odds Ratio", quotes = FALSE)
 odds             <- epitools::oddsratio(x=tab, method="wald") %>%
        print()
 or               <- odds$measure
clinical_sign     <- c(y,y)
  
  # if the assumptions of the chi squared test are violated...
  if (((cells_gte_five / total_cells)  < 0.8) || any(expected_counts < 1)) {
    
    fet           <- fisher.test(tab)
 

   p         <- fet[["p.value"]]
   
   selected <- ifelse(fet[["p.value"]]<=alpha, "yes", "no")
 
      or  <- cbind(clinical_sign , or, p, selected)%>%
      as.data.frame()%>%
     print()
    
      print("Fischer's test used", quotes = FALSE)
  } else {
    
     p  <- cst[["p.value"]]
     selected <- ifelse(cst[["p.value"]]<=alpha, "yes", "no")
   
      or <- cbind(clinical_sign , or, p, selected )%>%
      as.data.frame()%>%
      print()
      
      print("chi square test used", quotes = FALSE)
      }


}

```

## Clinical signs

-   Pale mucosae

-   Fever/hypothermia

-   Peripheral lymphadenopathy

-   Low BCS

-   Reduced appetite/anorexia

-   Ascites

-   Splenomegaly/ Hepatomegaly

-   Presence of wounds

Clinical signs of uncomplicated 'tick fever' - pale mucous membranes, reduced
appetite/ anorexia, fever, lymphadenomegaly, splenomegaly and hepatomegaly.

### *A. platys*

```{r}

uni_abmucsmemb_ana <- test_uni_cs(df=ticks_clean, y="cs_abmucsmemb",x="tbp_anaplasma")
uni_abtemp_ana <- test_uni_cs(df=ticks_clean, y="cs_abtemp", x="tbp_anaplasma")
uni_lowbcs_ana <- test_uni_cs(df=ticks_clean, y="cs_lowbcs", x="tbp_anaplasma")
uni_appetite_ana <- test_uni_cs(df=ticks_clean, y="cs_appetite", x="tbp_anaplasma")
uni_ascites_ana <- test_uni_cs(df=ticks_clean, y="cs_ascites", x="tbp_anaplasma")
uni_sphep_ana <- test_uni_cs(df=ticks_clean, y="cs_sphep", x="tbp_anaplasma")
uni_wound_ana <- test_uni_cs(df=ticks_clean, y="sk_wound", x="tbp_anaplasma")
uni_tickfever_ana <- test_uni_cs(df=ticks_clean, y="cs_tickfever", x="tbp_anaplasma")

```

### *E. canis*

```{r}

uni_abmucsmemb_ec <- test_uni_cs(df=ticks_clean, y="cs_abmucsmemb", x="tbp_ecanis")
uni_abtemp_ec <- test_uni_cs(df=ticks_clean, y="cs_abtemp", x="tbp_ecanis")
uni_lowbcs_ec <- test_uni_cs(df=ticks_clean, y="cs_lowbcs", x="tbp_ecanis")
uni_appetite_ec <- test_uni_cs(df=ticks_clean, y="cs_appetite", x="tbp_ecanis")
uni_ascites_ec <- test_uni_cs(df=ticks_clean, y="cs_ascites", x="tbp_ecanis")
uni_sphep_ec <- test_uni_cs(df=ticks_clean, y="cs_sphep", x="tbp_ecanis")
uni_wound_ec <- test_uni_cs(df=ticks_clean, y="sk_wound", x="tbp_ecanis")
uni_tickfever_ec <- test_uni_cs(df=ticks_clean, y="cs_tickfever",x="tbp_ecanis")

```

### Hemotropic mycoplasma spp.

```{r}

uni_abmucsmemb_myc <- test_uni_cs(df=ticks_clean, y="cs_abmucsmemb", x="tbp_mycoplasma")
uni_abtemp_myc <- test_uni_cs(df=ticks_clean, y="cs_abtemp", x="tbp_mycoplasma")
uni_lowbcs_myc <- test_uni_cs(df=ticks_clean, y="cs_lowbcs", x="tbp_mycoplasma")
uni_appetite_myc <- test_uni_cs(df=ticks_clean, y="cs_appetite", x="tbp_mycoplasma")
uni_ascites_myc <- test_uni_cs(df=ticks_clean, y="cs_ascites", x="tbp_mycoplasma")
uni_sphep_myc <- test_uni_cs(df=ticks_clean, y="cs_sphep", x="tbp_mycoplasma")
uni_wound_myc <- test_uni_cs(df=ticks_clean, y="sk_wound", x="tbp_mycoplasma")
uni_tickfever_myc <- test_uni_cs(df=ticks_clean, y="cs_tickfever", x="tbp_mycoplasma")

```

### *B. gibsoni*

```{r}

uni_abmucsmemb_bg <- test_uni_cs(df=ticks_clean, y="cs_abmucsmemb", x="tbp_bgibsoni")
uni_abtemp_bg <- test_uni_cs(df=ticks_clean, y="cs_abtemp", x="tbp_bgibsoni")
uni_lowbcs_bg <- test_uni_cs(df=ticks_clean, y="cs_lowbcs", x="tbp_bgibsoni")
uni_appetite_bg <- test_uni_cs(df=ticks_clean, y="cs_appetite", x="tbp_bgibsoni")
uni_ascites_bg <- test_uni_cs(df=ticks_clean, y="cs_ascites", x="tbp_bgibsoni")
uni_sphep_bg <- test_uni_cs(df=ticks_clean, y="cs_sphep", x="tbp_bgibsoni")
uni_wound_bg <- test_uni_cs(df=ticks_clean, y="sk_wound", x="tbp_bgibsoni")
uni_tickfever_bg <- test_uni_cs(df=ticks_clean, y="cs_tickfever", x="tbp_bgibsoni")

```

### *B. vogeli*

```{r}

uni_abmucsmemb_bv <- test_uni_cs(df=ticks_clean, y="cs_abmucsmemb", x="tbp_bvogeli")
uni_abtemp_bv <- test_uni_cs(df=ticks_clean, y="cs_abtemp", x="tbp_bvogeli")
uni_lowbcs_bv <- test_uni_cs(df=ticks_clean, y="cs_lowbcs", x="tbp_bvogeli")
uni_appetite_bv <- test_uni_cs(df=ticks_clean, y="cs_appetite", x="tbp_bvogeli")
uni_ascites_bv <- test_uni_cs(df=ticks_clean, y="cs_ascites", x="tbp_bvogeli")
uni_sphep_bv <- test_uni_cs(df=ticks_clean, y="cs_sphep", x="tbp_bvogeli")
uni_wound_bv <- test_uni_cs(df=ticks_clean, y="sk_wound", x="tbp_bvogeli")
uni_tickfever_bv <- test_uni_cs(df=ticks_clean, y="cs_tickfever", x="tbp_bvogeli")

```

### *H. canis*

```{r}

uni_abmucsmemb_hc <- test_uni_cs(df=ticks_clean, y="cs_abmucsmemb", x="tbp_hcanis")
uni_abtemp_hc <- test_uni_cs(df=ticks_clean, y="cs_abtemp", x="tbp_hcanis")
uni_lowbcs_hc <- test_uni_cs(df=ticks_clean, y="cs_lowbcs", x="tbp_hcanis")
uni_appetite_hc <- test_uni_cs(df=ticks_clean, y="cs_appetite", x="tbp_hcanis")
uni_ascites_hc <- test_uni_cs(df=ticks_clean, y="cs_ascites", x="tbp_hcanis")
uni_sphep_hc <- test_uni_cs(df=ticks_clean, y="cs_sphep", x="tbp_hcanis")
uni_wound_hc <- test_uni_cs(df=ticks_clean, y="sk_wound", x="tbp_hcanis")
uni_tickfever_hc <- test_uni_cs(df=ticks_clean, y="cs_tickfever", x="tbp_hcanis")

```

### Overall TBP infection

```{r}

uni_abmucsmemb_tbp <- test_uni_cs(df=ticks_clean, y="cs_abmucsmemb", x="tbp_infstat")
uni_abtemp_tbp <- test_uni_cs(df=ticks_clean, y="cs_abtemp", x="tbp_infstat")
uni_lowbcs_tbp <- test_uni_cs(df=ticks_clean, y="cs_lowbcs", x="tbp_infstat")
uni_appetite_tbp <- test_uni_cs(df=ticks_clean, y="cs_appetite", x="tbp_infstat")
uni_ascites_tbp <- test_uni_cs(df=ticks_clean, y="cs_ascites", x="tbp_infstat")
uni_sphep_tbp <- test_uni_cs(df=ticks_clean, y="cs_sphep", x="tbp_infstat")
uni_wound_tbp <- test_uni_cs(df=ticks_clean, y="sk_wound", x="tbp_infstat")
uni_tickfever_tbp <- test_uni_cs(df=ticks_clean, y="cs_tickfever", x="tbp_infstat")

```

### Mixed infections

```{r}

uni_abmucsmemb_infecnum <- test_uni_cs(df=ticks_clean, y="cs_abmucsmemb", x="tbp_infecnumgrp")
uni_abtemp_infecnum <- test_uni_cs(df=ticks_clean, y="cs_abtemp", x="tbp_infecnumgrp")
uni_lowbcs_infecnum <- test_uni_cs(df=ticks_clean, y="cs_lowbcs", x="tbp_infecnumgrp")
uni_appetite_infecnum <- test_uni_cs(df=ticks_clean, y="cs_appetite", x="tbp_infecnumgrp")
uni_ascites_infecnum <- test_uni_cs(df=ticks_clean, y="cs_ascites", x="tbp_infecnumgrp")
uni_sphep_infecnum <- test_uni_cs(df=ticks_clean, y="cs_sphep", x="tbp_infecnumgrp")
uni_wound_infecnum <- test_uni_cs(df=ticks_clean, y="sk_wound", x="tbp_infecnumgrp")
uni_tickfever_infecnum <- test_uni_cs(df=ticks_clean, y="cs_tickfever", x="tbp_infecnumgrp")

```

## Multivariable logistic regression models

#### Variables:

-   Age

-   Breed

-   Sex

-   Neuter status,

-   Tick infestation

-   Flea infestation

-   Geo-climatic zone

-   Use of ectoparasiticides

#### Saturated model :

[[pathogen infection] \~ age + breed + sex +neuter status+presence of
ticks+presence of fleas + geoclimatic zone + ectoparasiticide usage +
(1\|clinic)]{.smallcaps}

### *A. platys*

From univariable analysis:

```{r}

#significant variables

fit_glmm_anaplasma <- ticks_clean %>%
  lme4::glmer(
    tbp_anaplasma ~
        sk_tick+ 
    # sk_flea+
      dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_anaplasma)

tab_model(fit_glmm_anaplasma,
          show.intercept = TRUE,
          show.stat =  TRUE,
          show.se = TRUE,
          show.aic = TRUE)



#AIC start= 115.9
#AIC end= 115.7
```

All variables:

```{r}

# all variables
fit_glmm_anaplasma2 <- ticks_clean %>%
  lme4::glmer(
    tbp_anaplasma ~
   # dog_breedgrp2 +
   #dog_age_yr +
  #  dog_sex +
  # dog_neuterstat +
      sk_tick+ 
     #sk_flea+
   #   dog_tkrx +
    dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_anaplasma2)

tab_model(fit_glmm_anaplasma2)


#AIC start= 116.3
#AIC end= 115.7



```

### *E. canis*

From univariable analysis:

```{r}

#significant variables

fit_glmm_ecanis <- ticks_clean %>%
  lme4::glmer(
    tbp_ecanis ~
      #dog_breedgrp2 +
  # dog_age_yr +
   #dog_neuterstat +
      # sk_flea+
    # dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_ecanis)

tab_model(fit_glmm_ecanis,
          show.intercept = TRUE,
          show.est = TRUE,
          show.se = TRUE, 
          show.aic = TRUE)



#AIC start= 135.985
#AIC end=  
```

All variables:

```{r}

# all variables
fit_glmm_ecanis2 <- ticks_clean %>%
  lme4::glmer(
    tbp_ecanis ~
   # dog_breedgrp2 +
#   dog_age_yr +
    #dog_sex +
  # dog_neuterstat +
     # sk_tick+ 
     #sk_flea+
     # dog_tkrx +
   # dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_ecanis2)

tab_model(fit_glmm_ecanis2)


#AIC start= 116.3
#AIC end= 115.7



```

### Hemotropic mycoplasma spp.

From univariable analysis:

```{r}

#significant variables

fit_glmm_mycoplasma <- ticks_clean %>%
  lme4::glmer(
    tbp_mycoplasma ~
       dog_breedgrp2 +
     dog_sex +
        sk_flea+
     # dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_mycoplasma)

tab_model(fit_glmm_mycoplasma,
          show.stat = TRUE,
          show.intercept = TRUE,
          show.est = TRUE,
          show.se = TRUE,
          show.aic = TRUE)



#AIC start= 229.794
#AIC end= 230.659
```

All variables:

```{r}

# all variables
fit_glmm_mycoplasma2 <- ticks_clean %>%
  lme4::glmer(
    tbp_mycoplasma ~
    dog_breedgrp2 +
 #  dog_age_yr +
    dog_sex +
  # dog_neuterstat +
     # sk_tick+ 
     sk_flea+
     # dog_tkrx +
    #dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_mycoplasma2)

tab_model(fit_glmm_mycoplasma2,
          show.stat = TRUE,
          show.intercept = TRUE,
          show.est = TRUE,
          show.se = TRUE,
          show.aic = TRUE)



#AIC start=  195.1 
#AIC end= 230.659



```

### *B. gibsoni*

From univariable analysis:

```{r}

#significant variables

fit_glmm_bgibsoni <- ticks_clean %>%
  lme4::glmer(
    tbp_bgibsoni ~
dog_age_yr +
     #  dog_breedgrp2 +
     #  dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_bgibsoni)

tab_model(fit_glmm_bgibsoni,
          show.stat = TRUE,
          show.intercept = TRUE,
          show.est = TRUE,
          show.se = TRUE,
show.aic = TRUE)



#AIC start= 475.871
#AIC end= 474.749
```

All variables:

```{r}

# all variables
fit_glmm_bgibsoni2 <- ticks_clean %>%
  lme4::glmer(
    tbp_bgibsoni ~
   # dog_breedgrp2 +
   dog_age_yr +
   # dog_sex +
  # dog_neuterstat +
    #  sk_tick+ 
     #sk_flea+
     # dog_tkrx +
    #dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_bgibsoni2)

tab_model(fit_glmm_bgibsoni2,
show.aic = TRUE)


#AIC start= 381.708
#AIC end= 474.749



```

### *B. vogeli*

From univariable analysis:

```{r}

#significant variables

fit_glmm_bvogeli <- ticks_clean %>%
  lme4::glmer(
    tbp_bvogeli ~
#dog_sex +
      #sk_tick +
     #  dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_bvogeli)

tab_model(fit_glmm_bvogeli,
          show.intercept = TRUE,
          show.est = TRUE,
          show.se = TRUE,
show.aic = TRUE)



#AIC start= 140.098
#AIC end= 
```

All variables:

```{r}

# all variables
fit_glmm_bvogeli2 <- ticks_clean %>%
  lme4::glmer(
    tbp_bvogeli ~
   # dog_breedgrp2 +
   dog_age_yr +
   # dog_sex +
  # dog_neuterstat +
     # sk_tick+ 
     #sk_flea+
     # dog_tkrx +
   # dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_bvogeli2)

tab_model(fit_glmm_bvogeli2,
show.aic = TRUE)


#AIC start= 	120.800
#AIC end= 144.355



```

### *H. canis*

From univariable analysis:

```{r}

#significant variables

fit_glmm_hcanis <- ticks_clean %>%
  lme4::glmer(
    tbp_hcanis ~
dog_age_yr +
       sk_tick +
#  dog_neuterstat+
       dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_hcanis)

tab_model(fit_glmm_hcanis,
          show.stat = TRUE,
          show.intercept = TRUE,
          show.est = TRUE,
          show.se = TRUE,
show.aic = TRUE)



#AIC start= 392.135, 385.447
#AIC end= 406.233, 395.933
```

All variables:

```{r}

# all variables
fit_glmm_hcanis2 <- ticks_clean %>%
  lme4::glmer(
    tbp_hcanis ~
    #dog_breedgrp2 +
   dog_age_yr +
   # dog_sex +
  # dog_neuterstat +
      sk_tick+ 
   #  sk_flea+
     # dog_tkrx +
    dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_hcanis2)

tab_model(fit_glmm_hcanis2,
show.aic = TRUE)


#AIC start= 318.989
#AIC end= 	395.933



```

### Overall TBP infection

From univariable analysis:

```{r}

#significant variables

fit_glmm_infstat <- ticks_clean %>%
  lme4::glmer(
    tbp_infstat ~
  #dog_breedgrp2 +
   dog_age_yr +
    #dog_sex +
       # sk_tick+ 
   # dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_infstat)

tab_model(fit_glmm_infstat,
          show.intercept = TRUE,
          show.stat = TRUE,
          show.est = TRUE,
          show.se = TRUE,
show.aic = TRUE)



#AIC start= 449.614
#AIC end= 473.948
```

All variables:

```{r}

# all variables
fit_glmm_infstat2 <- ticks_clean %>%
  lme4::glmer(
    tbp_infstat ~
    #dog_breedgrp2 +
   dog_age_yr +
   # dog_sex +
   #dog_neuterstat +
     # sk_tick+ 
    # sk_flea+
     # dog_tkrx +
    #dog_zone2+
         (1|dog_clinic),
  data = .,
  family = "binomial"
  )
summary(fit_glmm_infstat2)

tab_model(fit_glmm_infstat2,
show.aic = TRUE)


#AIC start= 402.427	
#AIC end= 473.948


```

### Multivariable ordinal regression model for infection type

From univariable analysis:

```{r}

fit_clmm_infnum <- ordinal::clmm(tbp_infecnumgrp ~
                                # dog_breedgrp2 +
                           dog_age_yr +
                            dog_sex +
                            # sk_tick+ 
                              # dog_zone2+
                             (1|dog_clinic),
                         data = ticks_clean, link = "probit")

summary(fit_clmm_infnum)

tab_model(fit_clmm_infnum,
          show.intercept = TRUE,
          show.stat = TRUE,
          show.est = TRUE,
          show.se = TRUE,
show.aic = TRUE)


```

All variables:

```{r}

fit_clmm_infnum2 <- ordinal::clmm(tbp_infecnumgrp ~
                                 dog_breedgrp2 +
                           dog_age_yr +
                            dog_sex +
                           dog_neuterstat +
                           sk_tick+ 
                           sk_flea+
                           dog_tkrx +
                           dog_zone2+
                             (1|dog_clinic),
                         data = ticks_clean, link = "probit")

summary(fit_clmm_infnum2)

tab_model(fit_clmm_infnum2,
          show.intercept = TRUE,
          show.stat = TRUE,
          show.est = TRUE,
          show.se = TRUE,
show.aic = TRUE)


```
