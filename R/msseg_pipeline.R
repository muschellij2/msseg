#' @title MSSeg Pipeline
#' @description Full Pipeline for MSSeg Challenge
#' @param t1_pre T1-weighted pre-gad image/filename
#' @param t1_post T1-weighted post-gad image/filename
#' @param flair FLAIR image/filename
#' @param t2 T2-weighted image/filename
#' @param pd Proton Density pre-gad image/filename
#' @param gold_standard Gold Standar Lesion
#' @param outdir Output directory
#' @param outfile Output file for segmentation
#' @param num_templates Number of templates used in MASS
#' @param verbose print diagnostic messages
#' @param runcc Run connected components
#' @param min_vol Minimum volume of a lesion
#' @param remove_files Remove predictor files from hard disk
#' @param ... arguments passed to \code{\link{process_images}}
#' @return Final Segmentation
#' @export
#' @importFrom randomForest combine
#' @importFrom stats predict model.matrix
#' @importFrom ANTsR antsImageRead labelClusters as.array
#' @importFrom oro.nifti voxres
msseg_pipeline =  function(
  t1_pre = "3DT1.nii.gz",
  t1_post = "3DT1GADO.nii.gz",
  flair = "3DFLAIR.nii.gz",
  t2 = "T2.nii.gz",
  pd = "DP.nii.gz",
  gold_standard = NULL,
  outdir = ".",
  outfile = NULL,
  num_templates = 15,
  verbose = TRUE,
  runcc = TRUE,
  min_vol = 0.01,
  remove_files = FALSE,
  ...){

  process_images(t1_pre = t1_pre,
                 t1_post = t1_post,
                 flair = flair,
                 t2 = t2,
                 pd = pd,
                 gold_standard = gold_standard,
                 outdir = outdir,
                 num_templates = num_templates,
                 verbose = verbose, ...)
  # rm(list = "proc_L")
  # L = list(mask = ret_mask_fname,
  #          gold_standard = les_out_fname,
  #          all_imgs =   every_fname,
  #          outdir = outdir)

  if (verbose) {
    message("Making Processed Filenames")
  }
  have_gold_standard = !is.null(gold_standard)
  L = make_processed_filenames(
    types = c("FLAIR", "T1_Pre", "T1_Post", "T2", "PD"),
    gold_standard = have_gold_standard,
    outdir = outdir)

  if (verbose) {
    message("Making DF")
  }
  df = make_df(fnames = L$fnames,
               mask = L$mask,
               y_fname = L$gold_standard,
               verbose = verbose)

  if (verbose) {
    message("Reading Mask")
  }
  mask = check_nifti(L$mask)

  ind = which(df$qFLAIR > 0.5)

  p = rep(0, length = nrow(df))
  if (verbose) {
    message("Combining Random Forests")
  }
  rf = randomForest::combine(
    msseg::rf_model_1,
    msseg::rf_model_2,
    msseg::rf_model_3,
    msseg::rf_model_4,
    msseg::rf_model_5,
    msseg::rf_model_6,
    msseg::rf_model_7,
    msseg::rf_model_8,
    msseg::rf_model_9,
    msseg::rf_model_10)

  ##################################
  # Predict from the model
  ##################################
  if (verbose) {
    message("Predicting from Random Forest")
  }
  p[ind] = predict(rf,
                   newdata = df[ind,],
                   type = "prob")[, "1"]
  rm(list = "rf"); gc()
  rm(list = "df"); gc();
  if (verbose) {
    message("Remaking Image")
  }
  pimg = remake_img(p, img = mask, mask = mask)
  yhat = pimg > (msseg::cutoff)

  if (verbose) {
    message("Unstripping Image")
  }
  yhat = unstrip_image(strip_image = yhat,
                       outdir = outdir)
  if (!is.null(outfile)) {
    if (verbose) {
      message("Writing yhat")
    }
    writenii(yhat, filename = outfile)
  }
  if (runcc) {
    if (verbose) {
      message("Running Connected Components")
    }
    vres = voxres(yhat, units = "cm")
    aimg = antsImageRead(outfile)
    labs = labelClusters(aimg,
                         minClusterSize = 1,
                         fullyConnected = TRUE)
    min_vox = floor(min_vol / vres)
    tab = table(c(ANTsR::as.array(labs)))
    levs = names(tab[tab > min_vox])
    levs = as.numeric(levs)
    levs = levs[ levs > 0]
    olabs = ants2oro(labs)
    yhat = niftiarr(yhat, olabs %in% levs)
    if (!is.null(outfile)) {
      writenii(yhat, filename = outfile)
    }
  }
  if (remove_files) {
    fnames = L$fnames
    fnames = fnames[ file.exists(fnames) ]
    if (length(fnames) > 0) {
      file.remove(fnames)
    }
  }
  return(yhat)
}








