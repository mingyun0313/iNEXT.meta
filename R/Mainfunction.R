#' Data information for reference samples
#'
#' \code{DataInfobeta3Dmeta} provides basic data information in each combination of site and treatment for (1) the gamma reference sample in the pooled assemblage, and (2) the alpha reference sample in the joint assemblage.
#'
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a data.frame. The data frame has study/site and treatment as the first two columns, followed by columns for species names. Here an assemblage refers to a combination of study/site and treatment. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a data.frame. The data frame has study/site, treatment, patch as the first three columns, followed by columns for species names.
#' @param diversity Selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.
#' @param datatype Data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence data (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param PDtree (required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled data.
#' @param PDreftime (required only when \code{diversity = "PD"}), a numerical value specifying reference times for PD. Default is NULL (i.e., the age of the root of PDtree).
#' @param FDdistM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled data.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under a specified threshold value, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{FDtype = "AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_value"}), a numerical value between 0 and 1 specifying the tau value (threshold level) that will be used to compute FD. If \code{FDtau = "NULL"} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled data (i.e., quadratic entropy).
#'
#'
#'
#' @import iNEXT.beta3D
#' @import stringr
#'
#'
#' @return A data.frame including basic data information. \cr
#' For abundance data, basic information shared by TD, mean-PD and FD includes study/site name
#' (\code{Site}), treatment name (\code{Treatment}), pooled/joint assemblage (\code{Assemblage}), sample size (\code{n}), observed species
#' richness (\code{S.obs}), sample coverage estimates of the reference sample (\code{SC(n)}), sample coverage estimate
#' for twice the reference sample size (\code{SC(2n)}). Other additional information is given below. \cr
#'
#' \itemize{
#' (1) TD: the first five species abundance frequency counts in the reference sample (\code{f1–f5}).
#'
#' (2) Mean-PD: the the observed total branch length in the phylogenetic tree (\code{PD.obs}), the number of
#' singletons (\code{f1*}) and doubletons (\code{f2*}) in the node/branch abundance set, as well as the total branch
#' length of those singletons (\code{g1}) and of those doubletons (\code{g2}), and the reference time (\code{Reftime}).
#'
#' (3) FD (\code{FDtype = "AUC"}): the minimum distance (\code{dmin}) and the maximum distance (\code{dmax}) among
#' all non-diagonal elements in the distance matrix, and the mean distance between any two individuals
#' randomly selected from the dataset (\code{dmean}).
#'
#' (4) FD (\code{FDtype = "tau_value"}): the number of singletons (\code{a1*}) and of doubletons (\code{a2*}) among
#' the functionally indistinct set at the specified threshold level \code{'Tau'}, as well as the total contribution
#' of singletons (\code{h1}) and of doubletons (\code{h2}) at the specified threshold level \code{'Tau'}.
#' }
#'
#' For incidence data, the basic information for TD includes study/site name (\code{Site}), treatment name (\code{Treatment}),
#' pooled/joint assemblage (\code{Assemblage}), number of sampling units (\code{T}), total number of incidences (\code{U}), observed
#' species richness (\code{S.obs}), sample coverage estimates of the reference sample (\code{SC(T)}), sample coverage
#' estimate for twice the reference sample size (\code{SC(2T)}), as well as the first five species incidence
#' frequency counts (\code{Q1–Q5}) in the reference sample. For mean-PD and FD, output is similar to that
#' for abundance data.
#'
#'
#' @examples
#' ## (Data Information) Taxonomic diversity for abundance data
#' data("Spider_abundance_data")
#' info_TD_abu <- DataInfobeta3Dmeta(data = Spider_abundance_data, diversity = 'TD', datatype = 'abundance')
#' info_TD_abu
#'
#'
#' ## (Data Information) Taxonomic diversity for incidence data
#' data("Bat_incidence_data")
#' info_TD_inc <- DataInfobeta3Dmeta(data = Bat_incidence_data, diversity = 'TD', datatype = 'incidence_raw')
#' info_TD_inc
#'
#'
#' ## (Data Information) Mean phylogenetic diversity for abundance data
#' data("Spider_abundance_data")
#' data("Spider_tree")
#' info_PD_abu <- DataInfobeta3Dmeta(data = Spider_abundance_data, diversity = 'PD', datatype = 'abundance',
#'                                   PDtree = Spider_tree, PDreftime = NULL)
#' info_PD_abu
#'
#'
#' ## (Data Information) Mean phylogenetic diversity for incidence data
#' data("Bat_incidence_data")
#' data("Bat_tree")
#' info_PD_inc <- DataInfobeta3Dmeta(data = Bat_incidence_data, diversity = 'PD', datatype = 'incidence_raw',
#'                                   PDtree = Bat_tree, PDreftime = NULL)
#' info_PD_inc
#'
#'
#' ## (Data Information) Functional diversity for abundance data under a specified threshold level
#' data("Spider_abundance_data")
#' data("Spider_distM")
#' info_FDtau_abu <- DataInfobeta3Dmeta(data = Spider_abundance_data, diversity = 'FD', datatype = 'abundance',
#'                                      FDdistM = Spider_distM, FDtype = "tau_value", FDtau = NULL)
#' info_FDtau_abu
#'
#'
#' ## (Data Information) Functional diversity for abundance data when all threshold levels from 0 to 1 are
#' considered
#' data("Spider_abundance_data")
#' data("Spider_distM")
#' info_FDAUC_abu <- DataInfobeta3Dmeta(data = Spider_abundance_data, diversity = 'FD', datatype = 'abundance',
#'                                      FDdistM = Spider_distM, FDtype = "AUC", FDtau = NULL)
#' info_FDAUC_abu
#'
#'
#' ## (Data Information) Functional diversity for incidence data under a specified threshold level
#' data("Bat_incidence_data")
#' data("Bat_distM")
#' info_FDtau_inc <- DataInfobeta3Dmeta(data = Bat_incidence_data, diversity = 'FD', datatype = 'incidence_raw',
#'                                      FDdistM = Bat_distM, FDtype = "tau_value", FDtau = NULL)
#' info_FDtau_inc
#'
#'
#' ## (Data Information) Functional diversity for incidence data when all threshold levels from 0 to 1 are
#' considered
#' data("Bat_incidence_data")
#' data("Bat_distM")
#' info_FDAUC_inc <- DataInfobeta3Dmeta(data = Bat_incidence_data, diversity = 'FD', datatype = 'incidence_raw',
#'                                      FDdistM = Bat_distM, FDtype = "AUC", FDtau = NULL)
#' info_FDAUC_inc
#'
#' @export

DataInfobeta3Dmeta <- function(data, diversity = "TD", datatype = "abundance",
                               PDtree = NULL, PDreftime = NULL,
                               FDdistM = NULL, FDtype = "AUC", FDtau = NULL){

  if (datatype == "abundance"){

    data |>
      group_by(across(1:2)) |>
      summarise(mat = list(t(pick(everything())) %>%
                             DataInfobeta3D(data = ., diversity, datatype = "abundance", PDtree, PDreftime, FDdistM, FDtype, FDtau) |>
                             filter(Assemblage %in% c("Pooled assemblage", "Joint assemblage"))), .groups = "drop") |>
      unnest(mat) |>
      select(-Dataset)

  } else if (datatype == "incidence_raw") {

    data |>
      group_by(across(1:3)) |>
      summarise(mat = list(t(pick(everything()))), .groups = "drop") %>%
      group_by(across(1:2)) |>
      summarise(mat = list(mat), .groups = "keep") %>%
      reframe(mat %>%
                DataInfobeta3D(data = ., diversity, datatype = "incidence_raw", PDtree, PDreftime, FDdistM, FDtype, FDtau) |>
                filter(Assemblage %in% c("Pooled assemblage", "Joint assemblage"))) |>
      select(-Dataset)

  }

}


#' Estimates the difference of standardized 3D diversity with common sample coverage (for alpha, beta, gamma diversity, and four classes of dissimilarity measures) between two treatments, and fit a fixed- or random-effects model to perform meta analysis
#'
#' \code{iNEXTbeta3Dmeta} is a function that estimates the difference of standardized 3D (taxonomic, phylogenetic and functional) beta diversity between two treatments (e.g., enhanced vs. control), and perform meta analysis (fixed or random effect model) for several studies/sites.
#'
#' @param data (a) For \code{datatype = "abundance"}, data can be input as a data.frame. The data frame has study/site and treatment as the first two columns, followed by columns for species names. Here an assemblage refers to a combination of study/site and treatment. \cr
#' (b) For \code{datatype = "incidence_raw"}, data can be input as a data.frame. The data frame has study/site, treatment, patch as the first three columns, followed by columns for species names.
#' @param model Selection of model type: 'FE' = Fixed-effects model, 'RE' = Random-effects model.
#' @param diversity Selection of diversity type: 'TD' = Taxonomic diversity, 'PD' = Phylogenetic diversity, and 'FD' = Functional diversity.
#' @param order.q A numerical value specifying the diversity order, Default is \code{q = 0, 1, 2}.
#' @param datatype Data type of input data: individual-based abundance data (\code{datatype = "abundance"}) or species by sampling-units incidence data (\code{datatype = "incidence_raw"}) with all entries being 0 (non-detection) or 1 (detection).
#' @param level A numerical value between 0 and 1 specifying the sample coverage level used for computing standardized diversity and dissimilarity. By default(\code{level = NULL}), the function automatically calculates standardized 3D gamma, alpha, and beta diversities, along with four dissimilarity indices, up to the minimum coverage achieved by doubling the reference sample size across all site and treatment combinations.
#' @param nboot A positive integer specifying the number of bootstrap replications when assessing sampling uncertainty for estimating standardized beta3D diversity and the associated confidence intervals. Default is 10. If more accurate results are required, set \code{nboot = 100} (or \code{nboot = 200}).
#' @param treatment_order A character vector for the names of treatment. The difference of standardized 3D diversity will be computed as diversity of the first treatment minus the diversity of second treatment.
#' @param conf A positive number < 1 specifying the level of confidence interval. Default is 0.95.
#' @param PDtree (required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled data.
#' @param PDreftime (required only when \code{diversity = "PD"}), a numerical value specifying reference times for PD. Default is NULL (i.e., the age of the root of PDtree).
#' @param FDdistM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled data.
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_value"} for FD under a specified threshold value, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{FDtype = "AUC"}.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_value"}), a numerical value between 0 and 1 specifying the tau value (threshold level) that will be used to compute FD. If \code{FDtau = "NULL"} (default), then threshold is set to be the mean distance between any two individuals randomly selected from the pooled data (i.e., quadratic entropy).
#' @param FDcut_number (required only when \code{diversity = "FD"} and \code{FDtype = "AUC"}) a numeric number to cut [0, 1] interval into equal-spaced sub-intervals to obtain the AUC value by integrating the tau-profile. Equivalently, the number of tau values that will be considered to compute the integrated AUC value. Default is 30. A larger value can be set to obtain more accurate AUC value.
#'
#'
#'
#' @import iNEXT.beta3D
#' @import stringr
#'
#'
#' @return The function returns a list containing seven components: Gamma, Alpha, Beta, 1-C, 1-U, 1-V, and 1-S. Each component consists of two parts. The first part is a dataframe that includes a column named Site (or the user-defined name in the input data), representing the study or site identity. It also contains a column Difference, which indicates the difference in diversity between the two treatments,
#' along with SE for the standard error of the difference, LCL and UCL for the lower and upper confidence limits, and Order.q for the diversity order q. The type of diversity measure is recorded in the Diversity column, which can be TD, PD, or FD. Additionally, the dataframe includes two columns showing the estimated diversity values of the two treatments for each site,
#' and a Weight column representing the weight assigned to each site for the fixed or random effect model. The second part of each component is a summary table reporting meta-analytic statistics under the fixed or random effect model, including Cochran’s Q statistic (Q_val), the degrees of freedom (df_val), the associated p-value (p_val) for the heterogeneity test, the heterogeneity percentage (I2_val), and the estimated between-site variance (tau2_val).
#'
#' @examples
#'
#' ## Taxonomic diversity for abundance data
#'
#' # Coverage-based standardized TD
#' data("Spider_abundance_data")
#' output1_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "TD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
#'
#' output1_abu
#'
#' ## Phylogenetic diversity for abundance data
#'
#' # Coverage-based standardized PD
#' data("Spider_abundance_data")
#' data("Spider_tree")
#' output2_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "PD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                PDtree = Spider_tree, PDreftime = NULL, PDtype = "meanPD")
#' output2_abu
#'
#' ## Functional diversity for abundance data
#'
#' # Coverage-based standardized FD
#' data("Spider_abundance_data")
#' data("Spider_distM")
#' output3_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "FD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                FDdistM = Spider_distM, FDtype = "AUC", FDcut_number = 30)
#' output3_abu
#'
#' ## Taxonomic diversity for incidence data
#'
#' # Coverage-based standardized TD
#' data("Bat_incidence_data")
#' output1_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "TD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
#' output1_inc
#'
#' ## Phylogenetic diversity for incidence data
#'
#' # Coverage-based standardized PD
#' data("Bat_incidence_data")
#' data("Bat_tree")
#' output2_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "PD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                PDtree = Bat_tree, PDreftime = NULL, PDtype = "meanPD")
#' output2_inc
#'
#' ## Functional diversity for incidence data
#'
#' # Coverage-based standardized FD
#' data("Bat_incidence_data")
#' data("Bat_distM")
#' output3_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "FD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                FDdistM = Bat_distM, FDtype = "AUC", FDcut_number = 30)
#' output3_inc
#'
#' @export


iNEXTbeta3Dmeta <- function(data, model = "FE", diversity = "TD", order.q = 0, datatype = "abundance", level = NULL, nboot, treatment_order, conf = 0.95,
                            PDtree, PDreftime = NULL, PDtype = "meanPD",
                            FDdistM, FDtype = "AUC", FDtau = NULL, FDcut_number = 30){

  if (datatype == "abundance"){

    data[, 1] <- factor(data[, 1])
    name <- unique(data[, 1])

    if (is.null(level)){

      Cmin <- data |>
        group_by(across(c(1, 2))) |>
        group_split() |>
        lapply(function(x){

          x |>
            select(-c(1, 2)) |>
            unlist() %>%
            iNEXT.3D:::Coverage(., "abundance", 2 * sum(.))

        }) |>
        unlist() |>
        min()

    } else {

      Cmin <- level

    }

    data_T1 <- data %>%
      filter(.[[2]] == treatment_order[1]) %>%
      split(.[[1]]) |>
      lapply(function(x){
        x |> select(-c(1, 2)) |> t()
      })

    data_T2 <- data %>%
      filter(.[[2]] == treatment_order[2]) %>%
      split(.[[1]]) |>
      lapply(function(x){
        x |> select(-c(1, 2)) |> t()
      })

  } else if (datatype == "incidence_raw") {

    data[, 1] <- factor(data[, 1])
    name <- unique(data[, 1])

    if (is.null(level)){

      Cmin <- data |>
        group_by(across(1:3)) |>
        summarise(mat = list(t(pick(everything()))), .groups = "drop") |>
        group_by(across(1:2)) |>
        summarise(mat = list(do.call(rbind, mat)), .groups = "drop") %>%
        .$mat |>
        lapply(function(x) {

          iNEXT.3D:::Coverage(x, "incidence_raw", 2 * ncol(x))

        }) |>
        unlist() %>%
        min()
      # quantile(., 0.25)

    } else {

      Cmin <- level

    }

    data_T1 <- data %>%
      filter(.[[2]] == treatment_order[1]) %>%
      split(.[[1]]) %>%
      map(~ split(.x, .[[3]]) |>
            map(~ .x |> select(-c(1, 2, 3)) |> t())
      )

    data_T2 <- data %>%
      filter(.[[2]] == treatment_order[2]) %>%
      split(.[[1]]) %>%
      map(~ split(.x, .[[3]]) |>
            map(~ .x |> select(-c(1, 2, 3)) |> t())
      )

  } else {

    NULL

  }

  CC <- qnorm((1 - conf) / 2, lower.tail = F)

  ibeta_div_T1 <- iNEXTbeta3D(data_T1, diversity = diversity, q = order.q, datatype = datatype, level = Cmin, nboot = nboot,
                              PDtree = PDtree, PDreftime = PDreftime, PDtype = PDtype,
                              FDdistM = FDdistM, FDtype = FDtype, FDtau = FDtau, FDcut_number = FDcut_number)
  ibeta_div_T2 <- iNEXTbeta3D(data_T2, diversity = diversity, q = order.q, datatype = datatype, level = Cmin, nboot = nboot,
                              PDtree = PDtree, PDreftime = PDreftime, PDtype = PDtype,
                              FDdistM = FDdistM, FDtype = FDtype, FDtau = FDtau, FDcut_number = FDcut_number)

  gab_T1 <- ibeta_div_T1 |>
    map(~ imap(.x, ~ mutate(.x, Type = str_to_title(.y)) |>
                 rename_with(~ ifelse(. %in% c("Gamma", "Alpha", "Beta", "Dissimilarity"), "qD", .)))) |>
    map_dfr(~ bind_rows(.x))

  gab_T2 <- ibeta_div_T2 |>
    map(~ imap(.x, ~ mutate(.x, Type = str_to_title(.y)) |>
                 rename_with(~ ifelse(. %in% c("Gamma", "Alpha", "Beta", "Dissimilarity"), "qD", .)))) |>
    map_dfr(~ bind_rows(.x))

  finaldata <- lapply(c("Gamma", "Alpha", "Beta", "1-C", "1-U", "1-V", "1-S"), function(y) {

    div_T1 <- gab_T1 |> filter(Type == y)
    div_T2 <- gab_T2 |> filter(Type == y)

    merged_data_split <- merge(div_T1, div_T2, by = c("Dataset", "Order.q"), suffixes = c("_T1", "_T2")) %>%
      mutate(yi = qD_T1 - qD_T2, vi = s.e._T1^2 + s.e._T2^2) |>
      # escalc(measure = "MD",
      #        m1i = qD_T1, sd1i = s.e._T1, n1i = nboot,
      #        m2i = qD_T2,  sd2i = s.e._T2,  n2i = nboot,
      #        data = . |> mutate(nboot = nboot)) |>
      group_by(Order.q) |>
      mutate(Weight = ((1 / vi) / (sum(1 / vi))) * 100) |>
      ungroup() |>
      mutate(SE = sqrt(vi),
             LCL = yi - qnorm(1 - (1 - conf) / 2) * SE,
             UCL = yi + qnorm(1 - (1 - conf) / 2) * SE) |>
      rename(Difference = yi) |>
      group_split(Order.q)

    if (model == "FE"){

      results_list <- merged_data_split |>
        map(~ rma(yi = Difference, sei = SE, method = "FE", data = .x))

      meta_note <- map2_dfr(results_list, order.q, function(model, order_q) {
        tibble(
          Order.q = order_q,
          Q_val   = round(model$QE, 2),
          df_val  = model$k - 1,
          p_val   = ifelse(model$QEp < 0.001, "< 0.001", format(round(model$QEp, 3), nsmall = 3)),
          I2_val  = round(model$I2, 1),
          tau2_val= round(model$tau2, 3)
        )
      })

    } else if (model == "RE"){

      results_list <- merged_data_split |>
        map(~ rma(yi = Difference, sei = SE, method = "DL", data = .x))

      meta_note <- map2_dfr(results_list, order.q, function(model, order_q) {
        tibble(
          Order.q = order_q,
          Q_val   = round(model$QE, 2),
          df_val  = model$k - 1,
          p_val   = ifelse(model$QEp < 0.001, "< 0.001", format(round(model$QEp, 3), nsmall = 3)),
          I2_val  = round(model$I2, 1),
          tau2_val= round(model$tau2, 3)
        )
      })

    } else {

      stop('Please specify a valid model type: "FE" (Fixed Effects) or "RE" (Random Effects)')

    }

    meta_output <- lapply(1:length(order.q), function(x){

      results_list[[x]]$data |>
        select(Dataset, Difference, SE, LCL, UCL, Order.q, Diversity_T1, qD_T1, qD_T2, Weight) |>
        bind_rows(
          tibble(Dataset = paste(model, "Model"),
                 Difference = results_list[[x]]$beta |> as.numeric(),
                 SE = results_list[[x]]$se,
                 LCL = results_list[[x]]$ci.lb,
                 UCL = results_list[[x]]$ci.ub,
                 Order.q = results_list[[x]]$data$Order.q |> unique(),
                 Diversity_T1 = diversity,
                 qD_T1 = NA,
                 qD_T2 = NA,
                 Weight = 100)
        ) |>
        rename(Site = Dataset,
               Diversity = Diversity_T1) |>
        rename_with(~ treatment_order, .cols = c(qD_T1, qD_T2)) |>
        as.data.frame()

    }) |>
      map_df(~.x)

    return(list(meta_output, meta_note))

  })

  names(finaldata) <- paste("Summary", c("Gamma", "Alpha", "Beta", "1-C", "1-U", "1-V", "1-S"), sep = "_")
  return(finaldata)

}


#' Forest plot of the standardized 3D diversity difference between two treatments under fixed- or random-effects model
#'
#' \code{ggiNEXTmeta} is a function that provides forest plot for the difference of standardized 3D (taxonomic, phylogenetic and functional) beta diversity between two treatments.
#'
#' @param output The output of the iNEXTbeta3Dmeta function.
#' @param order.q A previously appeared 'Order.q' value in 'output'.
#' @param num_round A numerical value that the values show on the plot are rounded to the specified value of decimal places.
#' @param range Lower and upper limits for clipping confidence intervals to arrows.
#' @param type Specify diversity type (\code{'Gamma', 'Alpha', 'Beta'}), or dissimilarity type (\code{'1-C', '1-U', '1-V', '1-S'}).
#' @param level An optional sample coverage value (between 0 and 100 percent) to be annotated on the forest plot, indicating the fixed sample coverage used; if \code{level = NULL}, the annotation will be omitted.
#'
#' @import forestplot
#' @import grid
#'
#'
#' @return A forest plot that visualizing the output of iNEXTbeta3Dmeta. In the plot, it shows the difference of diversity between two treatments for each study/site and meta analysis (fixed- or random-effects model).
#'
#'
#' @examples
#' ## Taxonomic diversity for abundance data
#'
#' # Coverage-based standardized TD
#' data("Spider_abundance_data")
#' output1_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "TD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
#' ggiNEXTmeta(output1_abu, order.q = 0, num_round = 3, range = c(-20, 15), type = "Gamma", level = NULL)
#'
#'
#' ## Phylogenetic diversity for abundance data
#'
#' # Coverage-based standardized PD
#' data("Spider_abundance_data")
#' data("Spider_tree")
#' output2_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "PD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                PDtree = Spider_tree, PDreftime = NULL, PDtype = "meanPD")
#' ggiNEXTmeta(output2_abu, order.q = 0, num_round = 3, range = c(-20, 15), type = "Gamma", level = NULL)
#'
#'
#' ## Functional diversity for abundance data
#'
#' # Coverage-based standardized FD
#' data("Spider_abundance_data")
#' data("Spider_distM")
#' output3_abu <- iNEXTbeta3Dmeta(data = Spider_abundance_data, model = "RE", diversity = "FD",
#'                                order.q = c(0, 1, 2), datatype = "abundance", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                FDdistM = Spider_distM, FDtype = "AUC", FDcut_number = 30)
#' ggiNEXTmeta(output3_abu, order.q = 0, num_round = 3, range = c(-20, 15), type = "Gamma", level = NULL)
#'
#'
#' ## Taxonomic diversity for incidence data
#'
#' # Coverage-based standardized TD
#' data("Bat_incidence_data")
#' output1_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "TD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95)
#' ggiNEXTmeta(output1_inc, order.q = 0, num_round = 3, range = c(-20, 15), type = "Gamma", level = NULL)
#'
#'
#' ## Phylogenetic diversity for incidence data
#'
#' # Coverage-based standardized PD
#' data("Bat_incidence_data")
#' data("Bat_tree")
#' output2_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "PD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                PDtree = Bat_tree, PDreftime = NULL, PDtype = "meanPD")
#' ggiNEXTmeta(output2_inc, order.q = 0, num_round = 3, range = c(-20, 15), type = "Gamma", level = NULL)
#'
#'
#' ## Functional diversity for incidence data
#'
#' # Coverage-based standardized FD
#' data("Bat_incidence_data")
#' data("Bat_distM")
#' output3_inc <- iNEXTbeta3Dmeta(data = Bat_incidence_data, model = "RE", diversity = "FD",
#'                                order.q = c(0, 1, 2), datatype = "incidence_raw", level = NULL,
#'                                nboot = 10, treatment_order = c("Enhanced", "Control"), conf = 0.95,
#'                                FDdistM = Bat_distM, FDtype = "AUC", FDcut_number = 30)
#' ggiNEXTmeta(output3_inc, order.q = 0, num_round = 3, range = c(-20, 15), type = "Gamma", level = NULL)
#'
#'
#' @export


ggiNEXTmeta <- function(output, order.q, num_round = 3, range, type = NULL, level = NULL){

  if (!is.null(type)) output <- output[paste("Summary", type, sep = "_")]

  meta_output <- output[[1]][[1]]
  meta_note <- output[[1]][[2]]

  meta_output <- meta_output |> filter(Order.q == order.q)
  meta_note <- meta_note |> filter(Order.q == order.q)
  meta_note_q <- paste0(
    meta_output[nrow(meta_output), 1], " (Q = ", meta_note$Q_val,
    ", df = ", meta_note$df_val,
    ", p = ", meta_note$p_val,
    "; I² = ", meta_note$I2_val,
    ", τ² = ", meta_note$tau2_val, ")"
  )

  name <- colnames(meta_output)[1]
  name_treat <- colnames(meta_output)[c(8, 9)]
  ll <- length(meta_output[, 1])
  name_site <- meta_output[, 1][-ll]
  diversity <- meta_output$Diversity[1]
  order <- meta_output$Order.q[1]

  forestplot_q <- tibble(
    mean   = round(meta_output$Difference[-ll], num_round),
    lower  = round(meta_output$LCL[-ll], num_round),
    upper  = round(meta_output$UCL[-ll], num_round),
    study  = name_site,
    q_T1   = round(meta_output[, 8][-ll], num_round),
    q_T2   = round(meta_output[, 9][-ll], num_round),
    diff   = round(meta_output$Difference[-ll], num_round),
    LCL    = round(meta_output$LCL[-ll], num_round),
    UCL    = round(meta_output$UCL[-ll], num_round),
    w_fixed= paste(round(meta_output$Weight[-ll], 2), "%", sep = "")
  )

  fig <- forestplot_q |>
    forestplot(labeltext = c(study, q_T1, q_T2, diff, LCL, UCL, w_fixed),
               title = if (is.null(level)) {
                 type
               } else {
                 paste0(type, " (SC = ", level, "%)")
               },
               # clip = range,
               xlog = FALSE) |>
    fp_set_style(
      # box = "royalblue",
      # line = "darkblue",
      # summary = "royalblue",
      txt_gp = fpTxtGp(ticks = gpar(cex = 1))) |>
    fp_add_header(
      study = c("", name),
      q_T1  = c(paste(diversity, " (q = ", order, ")", sep = ""), name_treat[1]),
      q_T2  = c(paste(diversity, " (q = ", order, ")", sep = ""), name_treat[2]),
      diff  = c("", "Difference"),
      LCL   = c("", "LCL"),
      UCL   = c("", "UCL"),
      w_fixed = c("", "Weight")
    ) |>
    fp_append_row(
      mean  = round(meta_output$Difference[ll], num_round),
      lower = round(meta_output$LCL[ll], num_round),
      upper = round(meta_output$UCL[ll], num_round),
      study = "Meta analysis",
      diff  = round(meta_output$Difference[ll], num_round),
      LCL   = round(meta_output$LCL[ll], num_round),
      UCL   = round(meta_output$UCL[ll], num_round),
      is.summary = TRUE
    ) |>
    fp_decorate_graph(graph.pos = 4) |>
    fp_set_zebra_style("#EFEFEF") |>
    fp_add_lines()

  print(fig)
  grid.text(meta_note_q, x = 0.5, y = 0, just = "bottom", check.overlap = TRUE)

}

