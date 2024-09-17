#' eGFR by CKD-EPI equation
#'
#' A vectorised function to calculate estimated glomerular filtration rate using the CKD-EPI
#' equation. By default the equation accepts serum creatinine in µmol/l but can be changed to
#' mg/dl by setting the units parameter to "US". To allow for serial measurements over time, such as
#' for transplant follow-up data, there is an optional offset = n parameter which increases the age
#' value used in the equation by n years.
#'
#' Reference: Levey AS, Stevens LA, Schmid CH, et al. A new equation to estimate glomerular filtration
#' rate. Ann Intern Med 2009; 150(9):604-612.
#'
#' @param creat numeric vector of serum creatinine in µmol/l (or mg/dl if units = "US")
#' @param age numeric vector of age in years (accepts integers or decimals)
#' @param sex character vector of sex ("F" for female, "M" for male)
#' @param ethnicity character vector of patient ethnicity, one of "black" or "non-black"
#' @param units non-vectorised optional parameter for creatinine unit ("SI" for µmol/l (default), "US" for mg/dl)
#' @param offset non-vectorised optional numeric parameter for offset in years
#'
#' @return a numeric vector of eGFR values
#' @export
#'
#' @examples
#' ckd_epi(creat = 120, age = 45.2, sex = "M", ethnicity = "non-black")
#' ckd_epi(creat = 1.5, age = 64.3, sex = "F", ethnicity = "black", units = "US")
ckd_epi <- function(creat, age, sex, ethnicity, units = "SI", offset = 0) {

  # determine level of sexvar, alpha and kappa
  sexvar = ifelse(sex == "M", 1, 1.018)
  alpha = ifelse(sex == "M", -0.411, -0.329)
  kappa = ifelse(sex == "M", 0.9, 0.7)

  # convert SI to US units
  if (units == "SI") {
    creat = creat / 88.42
  }

  # age offset for serial measures
  if (offset > 0) age = age + offset

  gfr = 141 * ifelse(creat / kappa > 1, 1, creat / kappa)^alpha *
          ifelse(creat / kappa < 1, 1, creat / kappa)^-1.209 *
          0.993^age * sexvar

  # adjustment for black ethnicity
  gfr = ifelse(ethnicity == "black", gfr * 1.159, gfr)

  gfr
}

#' eGFR by abbreviated MDRD equation
#'
#' A vectorised function to calculate estimated glomerular filtration rate using the
#' abbreviated (four variable) MDRD equation. By default the equation accepts serum creatinine
#' in µmol/l but can be changed to mg/dl by setting the units parameter to "US".
#' To allow for serial measurements over time, such as
#' for transplant follow-up data, there is an optional offset = n parameter which increases the age
#' value used in the equation by n years.
#'
#' Reference: Levey AS, Greene T, Kusek JW, et al. A simplified equation to predict glomerular
#' filtration rate from serum creatinine. J Am Soc Nephrol 2000; 11:A0828.
#'
#' @param creat numeric vector of serum creatinine in µmol/l (or mg/dl if units = "US")
#' @param age numeric vector of age in years (accepts integers or decimals)
#' @param sex character vector of sex ("F" for female, "M" for male)
#' @param ethnicity character vector of patient ethnicity, one of "black" or "non-black"
#' @param units non-vectorised optional parameter for creatinine unit ("SI" for µmol/l (default), "US" for mg/dl)
#' @param offset non-vectorised optional parameter for offset in years
#'
#' @return a numeric vector of eGFR values
#' @export
#'
#' @examples
#' mdrd(creat = 120, age = 45.2, sex = "M", ethnicity = "non-black")
#' mdrd(creat = 1.5, age = 64.3, sex = "F", ethnicity = "black", units = "US")
mdrd <- function (creat, age, sex, ethnicity, units = "SI", offset = 0) {
  # determine level of sexvar
  sexvar = ifelse(sex == "M", 1, 0.742)

  # convert SI to US units
  if (units == "SI") {
    creat = creat / 88.42
  }

  # age offset for serial measures
  if (offset > 0) age = age + offset

  gfr =  186 * creat^-1.154 * age^-0.203 * sexvar

  # adjustment for black ethnicity
  gfr = ifelse(ethnicity == "black", gfr * 1.21, gfr)

  gfr
}

#' eGFR by bedside Schwartz formula
#'
#' A vectorised formula to calculate estimate glomerular filtration rate in children using the bedside Schwartz formula. By default this uses
#' serum creatinine in µmol/l but this can be changed to mg/dl by setting the optional units parameter to "US".
#'
#' Reference: Schwartz GJ, Munoz A, Schneider MF et al. New equations to estimate GFR in children
#' with CKD. J Am Soc Nephrol 2009; 20(3):629-637.
#'
#' @param creat numeric vector of creatinine levels in µmol/l (or mg/dl if units = "US")
#' @param height numeric vector of heights in cm
#' @param units non-vectorised optional parameter for creatinine unit ("SI" for µmol/l (default), "US" for mg/dl)
#'
#' @return numeric vector of eGFR values
#' @export
#'
#' @examples
#' # calculate using creatinine in µmol/l
#' schwartz(creat = 64, height = 101)
#'
#' # calculate using mg/dl
#' schwartz(creat = 0.7, height = 101, units = "US")
schwartz = function(creat, height, units = "SI") {
  if (units == "SI") {
    gfr = 36.5 * height / creat
  } else {
    gfr = 0.413 * height / creat
  }
  gfr
}

#' eGFR by CKD-EPI equation (US units)
#'
#' A wrapper function for the ckd_epi() vectorised function to calculate estimated glomerular
#' filtration rate using the CKD-EPI equation, using serum creatinine in mg/dl.
#' To allow for serial measurements over time, such as
#' for transplant follow-up data, there is an optional offset = n parameter which increases the age
#' value used in the equation by n years.
#'
#' Reference: Levey AS, Stevens LA, Schmid CH, et al. A new equation to estimate glomerular filtration
#' rate. Ann Intern Med 2009; 150(9):604-612.
#'
#' @param creat numeric vector of serum creatinine in µmol/l (or mg/dl if units = "US")
#' @param age numeric vector of age in years (accepts integers or decimals)
#' @param sex character vector of sex ("F" for female, "M" for male)
#' @param ethnicity character vector of patient ethnicity, one of "black" or "non-black"
#' @param offset non-vectorised optional parameter for offset in years
#'
#' @return a numeric vector of eGFR values
#' @export
#'
#' @examples
#' ckd_epi_US(creat = 1.5, age = 64.3, sex = "F", ethnicity = "black")
ckd_epi_US <- function (creat, age, sex, ethnicity, offset = 0) {
  ckd_epi(creat, age, sex, ethnicity, units = "US", offset = offset)
}

#' eGFR by abbreviated MDRD equation (US units)
#'
#' A wrapper for the mdrd4v() vectorised function to calculate estimated glomerular filtration rate
#' using the abbreviated (four variable) MDRD equation, but using serum creatinine in mg/dl.
#' To allow for serial measurements over time, such as
#' for transplant follow-up data, there is an optional offset = n parameter which increases the age
#' value used in the equation by n years.
#'
#' Reference: Levey AS, Greene T, Kusek JW, et al. A simplified equation to predict glomerular
#' filtration rate from serum creatinine. J Am Soc Nephrol 2000; 11:A0828.
#'
#' @param creat numeric vector of serum creatinine in µmol/l (or mg/dl if units = "US")
#' @param age numeric vector of age in years (accepts integers or decimals)
#' @param sex character vector of sex ("F" for female, "M" for male)
#' @param ethnicity character vector of patient ethnicity, one of "black" or "non-black"
#' @param offset non-vectorised optional parameter for offset in years
#'
#' @return a numeric vector of eGFR values
#' @export
#'
#' @examples

#' mdrd_US(creat = 1.5, age = 64.3, sex = "F", ethnicity = "black")
mdrd_US <- function (creat, age, sex, ethnicity, offset = 0) {
  mdrd(creat, age, sex, ethnicity, units = "US", offset = offset)
}

#' eGFR by bedside Schwartz formula (US units)
#'
#' A wrapper function for the schwartz() vectorised formula to calculate estimate glomerular filtration rate in children
#' using the bedside Schwartz formula, using
#' serum creatinine in mg/dl. Use the schwartz() function instead for µmol/l.
#'
#' Reference: Schwartz GJ, Munoz A, Schneider MF et al. New equations to estimate GFR in children
#' with CKD. J Am Soc Nephrol 2009; 20(3):629-637.
#'
#' @param creat numeric vector of creatinine levels in µmol/l (or mg/dl if units = "US")
#' @param height numeric vector of heights in cm
#'
#' @return numeric vector of eGFR values
#' @export
#'
#' @examples
#' # calculate using creatinine in -mg/dl
#' schwartz_US(creat = 0.7, height = 101)
schwartz_US = function(creat, height) {
  schwartz(creat = creat, height = height, units = "US")
}

#' Creatinine clearance by Cockcroft-Gault equation
#'
#' A vectorised function to estimate creatinine clearance using the Cockcroft-Gault equation.
#' By default this uses serum creatinine in µmol/l but can be changed to mg/dl by setting the
#' units parameter to "US"
#'
#' Reference: Cockcroft DW, Gault MH. Prediction of creatinine clearance from serum creatinine.
#' Nephron 1976; 16(1):31-41
#'
#' @param creat numeric vector of creatinine levels in µmol/l (or mg/dl if units = "US")
#' @param age numeric vector of ages in years
#' @param sex character vector of sex ("F" = female, "M" = male)
#' @param weight numeric vector of weights in kilograms
#' @param units non-vectorised parameter for creatinine units ("SI" for µmol/l (default) or "US" for mg/dl)
#'
#' @return numeric vector of creatinine clearances in ml/min
#' @export
#'
#' @examples
#' # calculate creatinine clearance using creatinine in µmol/l
#' cockcroft(creat = 88.4, age = 25, sex = "F", weight = 60)
#'
#' # calculate using creatinine in mg/dl
#' cockcroft(creat = 1, age = 25, sex = "F", weight = 60, units = "US")
cockcroft = function(creat, age, sex, weight, units = "SI"){
  if (units == "SI") {
    creat = creat / 88.4
  }

  sexvar = ifelse(sex == "M", 1, 0.85)
  sexvar * (140 - age) * weight / creat / 72
}


#' Creatinine clearance by Cockcroft-Gault equation (US units)
#'
#' A wrapper function for cockcroft(), a vectorised function to estimate creatinine clearance using the Cockcroft-Gault equation,
#' but using creatinine in mg/dl
#'
#' Reference: Cockcroft DW, Gault MH. Prediction of creatinine clearance from serum creatinine.
#' Nephron 1976; 16(1):31-41
#'
#' @param creat numeric vector of creatinine levels in mg/dl
#' @param age numeric vector of ages in years
#' @param sex character vector of sex ("F" = female, "M" = male)
#' @param weight numeric vector of weights in kilograms
#'
#' @return numeric vector of creatinine clearances in ml/min
#' @export
#'
#' @examples
#' cockcroft_US(creat = 1, age = 25, sex = "F", weight = 60)
cockcroft_US = function(creat, age, sex, weight){
  cockcroft(creat = creat, age = age, sex = sex, weight = weight, units = "US")
}

#' Ideal body weight
#'
#' A vectorised function to calculate adult ideal body weight based on height and sex.
#' This function assumes ideal BMI of 21.5 for females and 23 for males.
#'
#' @param height numeric vector of heights in cm
#' @param sex character vector of sex ("F" for female or "M" for male)
#'
#' @return numeric vector of ideal body weights in kg
#' @export
#'
#' @examples
#' ibw(height = 183, sex = "M")
ibw = function(height, sex) {
  bmivar = ifelse(sex == "M", 23, 21.5)
  height = height / 100
  height ^ 2 * bmivar
}


#' Creatinine based eGFR by CKiD U25 equation (US units)
#'
#' Reference: Pierce CB, Muñoz A, Ng DK, Warady BA, Furth SL, Schwartz GJ.
#' Age- and sex-dependent clinical equations to estimate glomerular filtration
#' rates in children and young adults with chronic kidney disease. Kidney
#' International. 2021;99(4):948–956. doi:10.1016/j.kint.2020.10.047
#'
#' @param creat numeric vector of creatinine levels in mg/dl
#' @param age numeric vector of ages in years
#' @param sex character vector of sex ("F" = female, "M" = male)
#' @param height numeric vector of heights in cm
#'
#' @return numeric vector of eGFR in ml/min/1.73m²
#' @export
#'
#' @examples
#' CKiD_U25_creatinine_US(creat = 1, age = 12, sex = "F", height = 132)
CKiD_U25_creatinine_US = function(creat, age, sex, height){

  if ( (min(age)<1) | (max(age)>25)) cat("\nWarning: there are age values <1 or >25 years; for those children, eGFR values might be invalid\n")

  coeff <- if(age < 12){
    if(sex == "F") 36.1 * 1.008 ^ (age - 12)
    else 39 * 1.008 ^ (age - 12)
  } else if (age < 18) {
    if(sex == "F") 36.1 * 1.023 ^ (age - 12)
    else 39 * 1.045 ^ (age - 12)
  } else {
    if(sex == "F") 41.4
    else 50.8
  }

  eGFR <- coeff * (height / 100) / creat
  round(eGFR, 1)
}

#' Cystatin C based eGFR by CKiD U25 equation (US units)
#'
#' Reference: Pierce CB, Muñoz A, Ng DK, Warady BA, Furth SL, Schwartz GJ.
#' Age- and sex-dependent clinical equations to estimate glomerular filtration
#' rates in children and young adults with chronic kidney disease. Kidney
#' International. 2021;99(4):948–956. doi:10.1016/j.kint.2020.10.047
#'
#' @param cystatin numeric vector of Cystatin C levels in mg/L
#' @param age numeric vector of ages in years
#' @param sex character vector of sex ("F" = female, "M" = male)
#'
#' @return numeric vector of eGFR in ml/min/1.73m²
#' @export
#'
#' @examples
#' CKiD_U25_cystatin_US(cystatin = 1, age = 18, sex = "F")
CKiD_U25_cystatin_US = function(cystatin, age, sex){

  if ( (min(age)<1) | (max(age)>25)) cat("\nWarning: there are age values <1 or >25 years; for those children, eGFR values might be invalid\n")

  coeff <- if(age < 12){
    if(sex == "F") 79.9*1.004 ^ (age - 12)
    else 87.2 * 1.011 ^ (age - 15)
  } else if(age < 15) {
    if(sex == "F") 79.9 * 0.974 ^ (age - 12)
    else 87.2 * 1.011 ^ (age-15)
  } else if (age < 18) {
    if(sex == "F") 79.9*0.974 ^ (age - 12)
    else 87.2 * 0.960 ^ (age - 15)
  } else {
    if(sex == "F") 77.1
    else 68.3
  }

  eGFR <- coeff / cystatin
  round(eGFR, 1)
}

#' Cystatin C and serum Creatinine based eGFR by CKiD U25 equation (US units)
#'
#' Wrapper function for estimating pediatric GFR based on data from the CKiD
#' study of children with CKD. Calculates a Cystatin C-based and Creatinine-
#' based eGFR then averages them. This has been shown to more faithfully
#' estimate measured GFR than either equation separately
#'
#' Reference: Pierce CB, Muñoz A, Ng DK, Warady BA, Furth SL, Schwartz GJ.
#' Age- and sex-dependent clinical equations to estimate glomerular filtration
#' rates in children and young adults with chronic kidney disease. Kidney
#' International. 2021;99(4):948–956. doi:10.1016/j.kint.2020.10.047
#'
#' @param cystatin numeric vector of Cystatin C levels in mg/dl
#' @param creat numeric vector of creatinine levels in mg/dl
#' @param age numeric vector of ages in years
#' @param sex character vector of sex ("F" = female, "M" = male)
#' @param height numeric vector of heights in cm
#' @param verbose a single boolean value. If true, return all the component
#' parts of the eGFR (cystatin, cr, combined), if false (defult), only return
#' the combined eGFR
#'
#' @return numeric vector of eGFR in ml/min/1.73m²
#' @export
#'
#' @examples
#' CKiD_U25_combined_US(cystatin = 1, creat = 0.7, age = 18, sex = "F",
#'                      height = 132, verbose = FALSE)
CKiD_U25_combined_US = function(cystatin, creat, age, sex, height,
                                verbose = FALSE){

  eGFRU25.cr <- CKiD_U25_creatinine_US(creat = creat, age = age,
                                    height = height, sex = sex)
  eGFRU25.cys <-  CKiD_U25_cystatin_US(cystatin = cystatin, age = age, sex = sex)
  sGFRU25.avg <- (eGFRU25.cys + eGFRU25.cr) / 2

  ## Determine if we should return all the component parts or just the result
  ## of the combined calculation
  if(verbose) {
    return(round(cbind(eGFRU25.cr, eGFRU25.cys, sGFRU25.avg),1))
  } else {
    return(round(sGFRU25.avg, 1))
  }
}
