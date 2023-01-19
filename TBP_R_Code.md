---
Title: "R code - TBP in Sri Lankan dogs"
author: "Ushani Atapattu"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float:
      collapsed: false
      smooth_scroll: false
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}

#Function to crosstabulate data 

tab_data <- function(df, x, y) { # df = data frame, x = variable 1, y = variable 2
  tab <-table(df[,x], df[,y], useNA = "ifany")
  print(tab)
}
```

## Summarising TBP infections

### TBP infection in the study group

```{r}
# Function to tabulate and calculate the confidence intervals

sum_data <- function(df, x, y){ # df = data frame, x = variable 1, y = variable 2
  
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

### Mixed TBP infections in the study group

```{r}

#tabulate mixed infection propotions

sum_data_type <- function(df, x, y){# df = data frame, x = variable 1, y = variable 2
  
# tabulate data in a contigency table
  tab <-table(df[ ,x], df[ ,y], useNA = "ifany")
  
  
#callculate the total from each group
  
 total <- marginSums(tab, margin = 1)
  
 
 tab <- cbind(tab, total)
 
 
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


## Univariable associations

### Univaruiable association between host, environment and other variable (except clinical signs)

```{r}
# Code for univariable association with fishers exact test or chi square test with significance at 0.2

test_uni     <-  function(df, x, y, alpha = 0.4) {
  
  #tabulate data in a contigency table
  
  tab             <- table(df[, x], df[, y])
  
 # print(tab)
  
  total_cells     <- nrow(tab) * ncol(tab)
  
  cst             <- chisq.test(tab, correct = FALSE)
  
  expected_counts <- cst$expected

  cells_gte_five  <- sum(expected_counts >= 5)
  
 # print("Odds Ratio", quotes = FALSE)
 odds             <- epitools::oddsratio(x=tab, method="wald")# %>%
        #print()
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
    
      #print("Fischer's test used", quotes = FALSE)
  } else {
    
     p  <- cst[["p.value"]]
     selected <- ifelse(cst[["p.value"]]<=alpha, "yes", "no")
   
      or <- cbind(variable, or, p, selected )%>%
      as.data.frame()%>%
      print()
      
      #print("chi square test used", quotes = FALSE)
      }


}

```


### Clinical Signs and TBP infection - Univariable Analysis

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


## Calculating p values, odds ratios, inverse logit and AIC values for ordered proportional odds logistic regression models


```{r}
model_values <- function(x){
  model <- summary(x) 
 
  
 # print(model)
coefficients <- model[["coefficients"]]
# print(coefficients)
  
 p_value <-   (1 - pnorm(abs(coefficients[ ,"t value"]), 0, 1))*2
  
odds_ratio <- exp(coefficients[ ,"Value"])
  
model_data <- cbind(coefficients, p_value, odds_ratio)
  
  ilogit <- exp(model[["zeta"]])/(1+exp(model[["zeta"]]))
  
print(noquote("Model values"))
  print(model_data)
  
  print(noquote("Inverse logit values"))
  
  print(ilogit)
  
   print(noquote("AIC"))
   
   print(AIC(x))
  
} 
