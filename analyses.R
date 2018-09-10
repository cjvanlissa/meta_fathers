# Converting effect sizes:
# https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf

poolSD <- function (ns, sds) { sqrt(sum((ns - 1) * sds)/(sum(ns) - length(ns))) }
g_calc <- function(m1, m2, n1, n2, sd1, sd2){ (m2-m1)/poolSD(c(n1,n2), c(sd1, sd2)) }
# Technically, the formula below applies the d to r conversion from 
# https://www.meta-analysis.com/downloads/Meta-analysis%20Converting%20among%20effect%20sizes.pdf
# however, Hedges' g is just d with a weighted pooled sd for unequal sample sizes
r_from_means <- function(m1, m2, n1, n2, sd1, sd2){
  a <- (n1+n2)^2/(n1*n2)
  d <- g_calc(m1, m2, n1, n2, sd1, sd2)
  d/sqrt(d^2+a)
}

set.seed(232)
data <- read.csv("Copy of 06-09-18_final_dataset2.csv", stringsAsFactors = FALSE)

# Drop empty columns added by Excel
rm_cols <- sapply(data, function(x){sum(is.na(x)) == nrow(data)})
data <- data[ , !rm_cols]

# Select only eligible cases
data <- data[data$sc_eligible == 1, ]

# Check which types of effect size exist
table(data$es_type)

# Convert mean differences to r
#data[data$es_type == "D", grep("es_", names(data))]
data[which(data$es_type == "D"), ]$es_correlation <- apply(data[which(data$es_type == "D"), c("es_mean1_lowinv", "es_mean2_highinv", "es_n1_lowinv", "es_n2_highinv", "es_s1_lowinv", "es_s2_highinv")], 1, function(x){ r_from_means(x[1], x[2], x[3], x[4], x[5], x[6]) })
data[which(data$es_type == "D"), ]$es_n <- apply(data[which(data$es_type == "D"), c("es_n1_lowinv", "es_n2_highinv")], 1, sum )
  
# Select cases with correlation effect size
analyzedat <- data[data$es_type %in% c("R", "Rs", "D", "R or F"), ]

# Check for missing effect size
analyzedat[which(is.na(analyzedat$es_correlation)), grep("es_", names(analyzedat))]
# Assume non-significant correlation is 0
analyzedat$es_correlation[which(is.na(analyzedat$es_correlation))] <- 0

# Check for missing effect size
sum(is.na(analyzedat$es_correlation))

# Drop redundant effect size variables
drop_es <- c("es_describe_result", "es_n1_lowinv", "es_n2_highinv", "es_mean1_lowinv", "es_mean2_highinv", "es_s1_lowinv", "es_s2_highinv", "es_beta", "es_b", "es_s_iv", "es_s_dv", "es_f_1df", "es_other")
analyzedat <- analyzedat[, -match(drop_es, names(analyzedat))]

# Drop unnecessary variables
drop_not_relevant <- c("coder_initials", 
                       #"sc_title",
                       "sc_journal", "sc_doi", "sc_scihub", "sc_corr_auth", "sc_email", "sc_eligible", "contact_authors", "sc_exclude_reason", "cor_1_with_2", "cor_1_with_3", "cor_1_with_4", "X")
analyzedat <- analyzedat[, -match(drop_not_relevant, names(analyzedat))]

# Drop redundant
drop_redundant <- c("sc_country")
analyzedat <- analyzedat[, -match(drop_redundant, names(analyzedat))]

# If reverse coded is missing, set it to 0 (we forgot to fill in some of these)
analyzedat$es_reverse_coded[is.na(analyzedat$es_reverse_coded)] <- 0
# If es_partial is missing, set it to 0 (we forgot to fill in some of these)
analyzedat$es_partial[is.na(analyzedat$es_partial)] <- 0
# If es_controls is missing, set it to 0 (we forgot to fill in some of these)
analyzedat$es_controls[is.na(analyzedat$es_controls)] <- 0

#table(analyzedat$sa_country, useNA = "always")

# Check which variables have the least missing
analyzedat[sapply(analyzedat, is.character)][analyzedat[sapply(analyzedat, is.character)] == ""] <- NA
analyzedat[sapply(analyzedat, is.character)][analyzedat[sapply(analyzedat, is.character)] == " "] <- NA
missingness <- apply(analyzedat, 2, function(x){sum(is.na(x))})
#71/377
#coder_initials sc_id sc_id_sample sc_id_effectsize sc_authors sc_pub_year sc_title sc_pub_type
#sc_eligible contact_authors sa_random sa_bias_demographics
#sa_bias_involvement sa_bias_language sa_age_m_iv sa_age_m_dv
table(analyzedat$sa_age_sd_dv)

analyzedat$sc_title[which(is.na(analyzedat$dv_timelag))]


drop_much_missing <- c("sa_larger_study")

# Note: Coley, Lewin-Bizan, & Carrano, the sc_id_sample and sc_id_effectsize might not be correctly coded

# Recode variables
tmp <- sapply(analyzedat$sa_data_year, gsub, pattern = "^(\\d{4})-\\d{4}$", replacement = "\\1")
tmp <- sapply(tmp, gsub, pattern = "^(\\d{4}).*?$", replacement = "\\1")
tmp[sapply(tmp, function(x){nchar(x)>4})] <- "1993"
analyzedat$sa_data_year <- as.integer(tmp)

numeric_vars <- c("sa_age_m_iv", "sa_age_sd_iv", "sa_age_m_dv", "sa_age_sd_dv")

proportion_vars <- c("sa_p_low_edu", "sa_p_low_income", "sa_p_low_ses", "sa_p_low_income_family")
analyzedat[proportion_vars] <- lapply(analyzedat[proportion_vars], function(x){
  x <- gsub("%", "", x)
  x[which(x > 1)] <- as.numeric(x[which(x > 1)])/100
  as.numeric(x)
})


# Impute missings

# For sa_p_low_ses_family, use highest available proportion SES indicator
analyzedat$sa_p_low_ses_family[which(is.na(analyzedat$sa_p_low_ses_family))] <- apply(
  analyzedat[which(is.na(analyzedat$sa_p_low_ses_family)), c("sa_p_low_edu", "sa_p_low_income", "sa_p_low_ses", "sa_p_low_income_family")],
  1, function(x){ ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)) })


# For missing data years, assume a publication lag of three years
analyzedat$sa_data_year[which(is.na(analyzedat$sa_data_year))] <- analyzedat$sc_pub_year[which(is.na(analyzedat$sa_data_year))]-3

# For missing IV reliability, try to find similar scales
analyzedat$iv_instrumentname[which(is.na(analyzedat$iv_reliability))]
#analyzedat$iv_instrumentname[which(is.na(analyzedat$iv_reliability))][c(1, 6, 26, 27, 33, 34)]

# Probably better NOT to set single item reliability to 1; instead, assign average reliability
#single_item <- c("no instrument, mother report",
#                 "‘How many days has father seen child during the past 30 days?",
#                 "(i)n the past month, have you spanked (child) because (he/she) was misbehaving or acting up?’ (Fragile Families, 2005). Parents who reported spanking were also asked about frequency of spanking.",
#                 "how often they read to their child on a 6-point Likert scale:rarely, not at all, a few times a month, a few times a week, once aday, or more than once a day.",
#                 "language activity time using 5-point scales (1 = rarely or never, 2 = occasionally, but not on a regular basis, 3 = 15-20 minutes per week, 4 = 15-20 minutes several times per week, 5 = 15-20 minutes each day)",
#                 "Interaction time in hours")
#analyzedat$iv_reliability[which(analyzedat$iv_instrumentname %in% single_item)] <- 1

analyzedat$iv_reliability[which(grepl("SALT", toupper(analyzedat$iv_instrumentname)) & is.na(analyzedat$iv_reliability))] <- 
  mean(as.numeric(analyzedat$iv_reliability[grepl("SALT", toupper(analyzedat$iv_instrumentname))]), na.rm = TRUE)

analyzedat$iv_reliability[which(grepl("CHILDES", toupper(analyzedat$iv_instrumentname)) & is.na(analyzedat$iv_reliability))] <- 
  mean(as.numeric(analyzedat$iv_reliability[grepl("CHILDES", toupper(analyzedat$iv_instrumentname))]), na.rm = TRUE)

# For missing DV reliability, try to find similar scales
#analyzedat$dv_instrumentname[which(is.na(analyzedat$dv_reliability))]
#any(is.na(analyzedat$dv_reliability))

analyzedat$dv_reliability[which(grepl("SALT", toupper(analyzedat$dv_instrumentname)) & is.na(analyzedat$dv_reliability))] <- 
  mean(as.numeric(analyzedat$dv_reliability[grepl("SALT", toupper(analyzedat$dv_instrumentname))]), na.rm = TRUE)

analyzedat$dv_reliability[which(grepl("CHILDES", toupper(analyzedat$dv_instrumentname)) & is.na(analyzedat$dv_reliability))] <- 
  mean(as.numeric(analyzedat$dv_reliability[grepl("CHILDES", toupper(analyzedat$dv_instrumentname))]), na.rm = TRUE)

# PLAI : .94, Test Review: Preschool Language Assessment Instrument (PLAI) https://sites.ualberta.ca/~lphillip/documents/Preschool%20Language%20Assessment%20Instrument-2%20(PLAI-2).doc
analyzedat$dv_reliability[which(grepl("PLAI", toupper(analyzedat$dv_instrumentname)) & is.na(analyzedat$dv_reliability))] <- 
.94

# Bowles, Pentimonti, Gerde, & Montroy, 2013: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.1009.1424&rep=rep1&type=pdf
analyzedat$dv_reliability[which(grepl("^IRT-BASED", toupper(analyzedat$dv_instrumentname)) & is.na(analyzedat$dv_reliability))] <- 
  .94

analyzedat$dv_reliability[which(grepl("SICD", toupper(analyzedat$dv_instrumentname)) & is.na(analyzedat$dv_reliability))] <- 
  mean(as.numeric(analyzedat$dv_reliability[grepl("SICD", toupper(analyzedat$dv_instrumentname))]), na.rm = TRUE)

analyzedat$dv_reliability[which(grepl("REYNELL", toupper(analyzedat$dv_instrumentname)) & is.na(analyzedat$dv_reliability))] <- 
  mean(as.numeric(analyzedat$dv_reliability[grepl("REYNELL", toupper(analyzedat$dv_instrumentname))]), na.rm = TRUE)

analyzedat$dv_reliability[which(grepl("(PEABODY|PPVT)", toupper(analyzedat$dv_instrumentname)) & is.na(analyzedat$dv_reliability))] <- 
  mean(as.numeric(analyzedat$dv_reliability[grepl("(PEABODY|PPVT)", toupper(analyzedat$dv_instrumentname))]), na.rm = TRUE)

analyzedat$dv_reliability[which(grepl("(MA?CARTHUR|CDI)", toupper(analyzedat$dv_instrumentname)) & is.na(analyzedat$dv_reliability))] <- 
  mean(as.numeric(analyzedat$dv_reliability[grepl("(MA?CARTHUR|CDI)", toupper(analyzedat$dv_instrumentname))]), na.rm = TRUE)




#analyzedat$dv_reliability[!grepl("^[0-9]*\\.?[0-9]+$", analyzedat$dv_reliability)]
# 
# table(analyzedat$sa_larger_study[which(is.na(analyzedat$sa_age_m_f))], useNA = "always")
# x <- analyzedat$sa_age_m_f
# table(x)
# sum(is.na(x))


# Categorize manually categorized variables

categorize <- c("iv_type", "dv_instrumentname")
for(var in categorize){
  categories <- read.csv(paste0("categorize_", var, ".csv"))
  analyzedat[[paste0(var, "_cat")]] <- categories[match(analyzedat[[var]], categories$Original), ]$Category
}


# Categorize manually categorized variables

categorize <- c("sa_larger_study", "iv_type", "dv_instrumentname")
for(var in categorize){
  categories <- read.csv(paste0("categorize_", var, ".csv"), stringsAsFactors = FALSE)
  analyzedat[[paste0(var, "_cat")]] <- categories[match(analyzedat[[var]], categories$Original), ]$Category
}

# Lowercase to uppercase
lowercase_vars <- c("iv_location", "sa_pop_dens")
analyzedat[lowercase_vars] <- lapply(lowercase_vars, function(x) toupper(analyzedat[[x]]))


# For missing sample variables, try to use other studies from the same sample
sample_vars_cat <- c("sa_random", "sa_bias_demographics", "sa_bias_involvement", "sa_bias_language", "sa_country", "sa_pop_dens", "sa_mother")
sample_vars_num <- c("sa_data_year", "sa_approached_participated", "sa_began_finished", "sa_age_m_iv", "sa_age_sd_iv", "sa_age_m_dv", "sa_age_sd_dv", "sa_p_male", "sa_age_m_f", "sa_age_sd_f", "sa_p_single_parent", "sa_p_hispanic_f", "sa_p_black_f", "sa_p_asian_f", "sa_p_otherethnic_f", "sa_p_low_ses_family", "sa_p_bilingual", "sa_p_intact", "sa_p_firstborn")
studies <- names(table(analyzedat$sa_larger_study_cat)[-1])

analyzedat$sa_larger_study_cat[is.na(analyzedat$sa_larger_study_cat)] <- paste(analyzedat$sc_authors[is.na(analyzedat$sa_larger_study_cat)], analyzedat$sc_pub_year[is.na(analyzedat$sa_larger_study_cat)], sep = ", ")

# Impute numeric
# Check missings
apply(analyzedat[sample_vars_num], 2, function(x){sum(is.na(x))})
for(var in sample_vars_num){
  for(study in studies){
    calc_average <- mean(as.numeric(analyzedat[[var]][analyzedat$sa_larger_study_cat == study & !is.na(analyzedat[[var]])]), na.rm = TRUE)
    analyzedat[[var]][analyzedat$sa_larger_study_cat == study & is.na(analyzedat[[var]])] <- ifelse(is.nan(calc_average), NA, calc_average)
  }
}
# Check missings
apply(analyzedat[sample_vars_num], 2, function(x){sum(is.na(x))})

# Impute categorical
# Check missings
apply(analyzedat[sample_vars_cat], 2, function(x){sum(is.na(x))})
for(var in sample_vars_cat){
  for(study in studies){
    tab <- table(analyzedat[[var]][analyzedat$sa_larger_study_cat == study & !is.na(analyzedat[[var]])])
    tab <- names(tab)[which.max(tab)]
    analyzedat[[var]][analyzedat$sa_larger_study_cat == study & is.na(analyzedat[[var]])] <- ifelse(is.null(tab), NA, tab)
  }
}
# Check missings
apply(analyzedat[sample_vars_cat], 2, function(x){sum(is.na(x))})


# Create new ID variables
analyzedat$sc_auth_year <- paste(analyzedat$sc_authors, analyzedat$sc_pub_year, sep = ", ")

analyzedat$sc_auth_year[which(analyzedat$sc_auth_year %in% unique(analyzedat$sc_auth_year[analyzedat$sc_id_sample == 2]))] <-
  paste(analyzedat$sc_auth_year[which(analyzedat$sc_auth_year %in% unique(analyzedat$sc_auth_year[analyzedat$sc_id_sample == 2]))], analyzedat$sc_id_sample[which(analyzedat$sc_auth_year %in% unique(analyzedat$sc_auth_year[analyzedat$sc_id_sample == 2]))], sep = ", Sample ")

analyzedat$id_study <- as.numeric(as.factor(analyzedat$sc_auth_year))
analyzedat <- analyzedat[order(analyzedat$id_study), ]

analyzedat$id_row <- 1:nrow(analyzedat)

# Remove variables no longer needed (e.g., only used for checks and computations)
drop_redundant <- c("sc_id", "sc_id_sample", "sc_id_effectsize", "sc_title", "sa_larger_study",
                    "sa_educationlevel_f", "sa_p_low_edu",
                  "sa_income_f", "sa_p_low_income", "sa_ses_f", "sa_p_low_ses",
                  "sa_income_family", "sa_p_low_income_family", "sa_ses_family", "iv_type",
                  "dv_des", "dv_instrumentname", "sc_auth_year",
                  "iv_des", "iv_instrumentname")
analyzedat <- analyzedat[, -match(drop_redundant, names(analyzedat))]

numeric_vars <- c("sc_impact_factor_that_year", "sa_data_year", "sa_approached_participated", "sa_began_finished", "sa_p_male", "sa_p_single_parent", "sa_p_hispanic_f", "sa_p_black_f", "sa_p_asian_f", "sa_p_otherethnic_f", "sa_p_low_ses_family", "sa_p_bilingual", "sa_p_bio", "sa_p_intact", "sa_p_firstborn", "es_controls", "es_correlation", "id_study", "iv_reliability", "dv_timelag", "dv_reliability", "sa_age_m_iv", "sa_age_sd_iv", "sa_age_m_dv", "sa_age_sd_dv", "sa_age_m_f", "sa_age_sd_f", "sc_pub_year", "es_n", "id_row")
categorical_vars <- c("es_reverse_coded", "es_partial", "sc_pub_type", "sa_random", "sa_bias_demographics", "sa_bias_involvement", "sa_bias_language", "sa_country", "sa_pop_dens", "sa_mother", "iv_assessment", "iv_location", "iv_informant", "dv_type", "dv_assessment", "dv_location", "dv_informant", "es_type", "iv_type_cat", "dv_instrumentname_cat", "sa_larger_study_cat", "de_experiment", "de_longitudinal", "de_manipulation", "de_larger_study", "iv_qualitative", "iv_researcher_present", "iv_recording", "iv_validity", "dv_researcher_present", "dv_recording", "dv_validity")

analyzedat[numeric_vars] <- lapply(analyzedat[numeric_vars], as.numeric)
analyzedat[categorical_vars] <- lapply(analyzedat[categorical_vars], as.factor)

library(missForest)

impute_data <- missForest(analyzedat[, -match(c("id_study", "id_row", "sc_authors"), names(analyzedat))], variablewise = TRUE)
impute_data$OOBerror

analyzedat_imputed <- data.frame(analyzedat[, match(c("id_study", "id_row", "sc_authors"), names(analyzedat))], impute_data$ximp)

write.csv(analyzedat_imputed, "imputed_data.csv", row.names = FALSE)