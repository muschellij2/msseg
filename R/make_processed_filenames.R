
#' Make Processed Filenames
#' @description Get filename output from processing
#' @param types Images used
#' @param gold_standard Should the gold standard be added
#' @param outdir output directory of processing
#'
#' @return data.frame of values
#' @export
make_processed_filenames = function(
  types = c("FLAIR", "T1_Pre", "T1_Post", "T2", "PD"),
  gold_standard = TRUE,
  outdir = "."){

  wide = data.frame(id = 1,
                    stringsAsFactors = FALSE)

  fnames = file.path(outdir, paste0(types, ".nii.gz"))
  fnames = t(fnames)
  colnames(fnames) = types
  df = data.frame(fnames,
                  stringsAsFactors = FALSE)
  wide = cbind(wide, df)
  rm(list = c("df", "fnames"))

  wide$id = NULL
  # wide$proc_dir = outdir
  if (gold_standard) {
    y_fname = file.path(outdir,
                        "GoldStandard_strip.nii.gz")
  } else {
    y_fname = NULL
  }
  wide$mask = file.path(outdir,
                        "FLAIR_MALF_Brain_Mask_strip.nii.gz")



  #######################################
  # Resetting the dir to processed
  ########################################
  gg = function(x, norm = TRUE,
                quant = FALSE){

    x = nii.stub(x, bn = TRUE)
    x = gsub("\\d\\d\\d_(.*)", "\\1", x)
    x = paste0(x, "_strip")
    if (norm) {
      x = paste0(x, "_norm")
    }
    if (quant) {
      x = paste0(x, "_quantile")
    }
    x = paste0(x, ".nii.gz")
    x = file.path(outdir, x)
    x
  }
  wide$qFLAIR = gg(wide$FLAIR,
                   norm = FALSE,
                   quant = TRUE)

  wide$qT1_Pre = gg(wide$T1_Pre,
                    norm = FALSE,
                    quant = TRUE)
  wide$qT1_Post = gg(wide$T1_Post,
                     norm = FALSE,
                     quant = TRUE)
  wide$qT2 = gg(wide$T2,
                norm = FALSE,
                quant = TRUE)
  wide$qPD = gg(wide$PD,
                norm = FALSE,
                quant = TRUE)


  wide$FLAIR = gg(wide$FLAIR)
  wide$T1_Pre = gg(wide$T1_Pre)
  wide$T1_Post = gg(wide$T1_Post)
  wide$T2 = gg(wide$T2)
  wide$PD = gg(wide$PD)

  # if (gold_standard) {
  #   wide$GoldStandard = file.path(
  #     outdir,
  #     "GoldStandard_strip.nii.gz")
  # }

  wide$flair15 = file.path(
    outdir,
    "Over_FLAIR_15Quantile.nii.gz")

  wide$flair90 = file.path(
    outdir,
    "FLAIR_90.nii.gz")

  keep_cols = types
  icol = keep_cols[1]

  for (icol in keep_cols) {
    all_fnames = wide[, icol]
    all_fname_stubs = nii.stub(all_fnames)

    ################################
    # Flipped Differences
    ################################
    fnames = paste0(all_fname_stubs,
                    "_Flipped_Difference.nii.gz")
    df = data.frame(fnames,
                    stringsAsFactors = FALSE)
    colnames(df) = paste0(icol, "_flip")
    wide = cbind(wide, df)

    ##############################
    # Smoothed Images with 3mm
    # Gaussian and write out
    ##############################
    sigmas = c(3, 10, 20)
    fnames = outer(all_fname_stubs,
                   paste0("_smooth", sigmas, ".nii.gz"),
                   paste0)
    df = data.frame(fnames,
                    stringsAsFactors = FALSE)
    colnames(df) = paste0(icol, "_smooth",
                          sigmas)
    wide = cbind(wide, df)

    ##############################
    # Perona Malik ones
    ##############################
    fnames = paste0(all_fname_stubs,
                    "_PeronaMalik.nii.gz")
    df = data.frame(fnames,
                    stringsAsFactors = FALSE)
    colnames(df) = paste0(icol,
                          "_PeronaMalik")
    wide = cbind(wide, df)

    ################################
    # Z to a template
    ################################
    zs = c("_Z_Median", "_Z_Mean")
    fnames = outer(all_fname_stubs,
                   paste0(zs, ".nii.gz"),
                   paste0)
    df = data.frame(fnames,
                    stringsAsFactors = FALSE)
    colnames(df) = paste0(icol, tolower(zs))
    wide = cbind(wide, df)

    ################################
    # Z to a template
    ################################
    zs = c("_Z_Median", "_Z_Mean")
    zs = paste0("_DRAMMS", zs)
    fnames = outer(all_fname_stubs,
                   paste0(zs, ".nii.gz"),
                   paste0)
    df = data.frame(fnames,
                    stringsAsFactors = FALSE)
    colnames(df) = paste0(icol, tolower(zs))
    wide = cbind(wide, df)

    ################################
    # Local Moments
    ################################
    moments =  c("mn", "sd", "skew",
                 "kurt", "grad", "z")
    fnames = outer(paste0(all_fname_stubs,
                          "_"),
                   paste0(moments, ".nii.gz"),
                   paste0)
    df = data.frame(fnames,
                    stringsAsFactors = FALSE)
    colnames(df) = paste0(icol, "_",
                          tolower(moments))
    wide = cbind(wide, df)
  }

  wide$class = file.path(outdir,
                         paste0("T1_Pre_MALF_Tissue_Classes",
                                ".nii.gz") )
  wide$gclass = file.path(outdir,
                          paste0("T1_Pre_MALF_Tissue_Classes",
                                 "_Gauss", ".nii.gz") )
  wide$ants_seg = file.path(outdir,
                            "T1_Pre_ANTs_Tissue_Classes.nii.gz")
  wide$flair_ants_seg = file.path(
    outdir,
    "FLAIR_ANTs_Tissue_Classes.nii.gz")

  tissue_classes = 1:3
  probs = outer(outdir,
                paste0("T1_Pre_ANTs_Tissue_Prob_",
                       tissue_classes, ".nii.gz"),
                file.path)
  colnames(probs) = paste0("prob_",
                           tissue_classes)
  df = data.frame(probs,
                  stringsAsFactors = FALSE)
  wide = cbind(wide, df)
  rm(list = c("df", "probs"))

  tissue_classes = 1:4
  probs = outer(outdir,
                paste0("FLAIR_ANTs_Tissue_Prob_",
                       tissue_classes, ".nii.gz"),
                file.path)
  colnames(probs) = paste0("flair_prob_",
                           tissue_classes)
  df = data.frame(probs,
                  stringsAsFactors = FALSE)
  wide = cbind(wide, df)
  rm(list = c("df", "probs"))

  probs = file.path(outdir,
                    paste0(types, "_strip_norm",
                           "_self_flip.nii.gz")  )
  probs = t(probs)
  colnames(probs) = paste0(types, "_self_flip")
  df = data.frame(probs,
                  stringsAsFactors = FALSE)

  wide = cbind(wide, df)
  rm(list = c("df", "probs"))

  all_names = c("T1_Pre", "FLAIR",
                "T1_Post",
                "T2",
                "PD")
  eg = expand.grid(img1 = "FLAIR",
                   img2 = setdiff(all_names, "FLAIR"))
  eg = as.matrix(eg)
  df = as.data.frame(eg,
                     stringsAsFactors = FALSE)
  df$iname = paste0(df$img1, "_", df$img2,
                    "_Correlation")
  df$name = paste0(df$iname,
                   ".nii.gz")
  fnames = outer(outdir,
                 df$name,
                 file.path)
  colnames(fnames) = tolower(df$iname)
  fnames = as.data.frame(fnames,
                         stringsAsFactors = FALSE)
  wide = cbind(wide, fnames)

  mask_fname = file.path(outdir,
                         "FLAIR_MALF_Brain_Mask_strip.nii.gz")
  L = list(mask = mask_fname,
           fnames = unlist(wide))
  L$gold_standard = y_fname
  return(L)
}

