#' @title MSSeg Pipeline
#' @description Full Pipeline for MSSeg Challenge
#' @param t1_pre T1-weighted pre-gad image/filename
#' @param t1_post T1-weighted post-gad image/filename
#' @param flair FLAIR image/filename
#' @param t2 T2-weighted image/filename
#' @param pd Proton Density pre-gad image/filename
#' @param gold_standard Gold Standar Lesion
#' @param outdir Output directory
#' @param num_templates Number of templates used in MASS
#' @param verbose print diagnostic messages
#' @return Final Segmentation
#' @export
msseg_pipeline =  function(
  t1_pre = "3DT1.nii.gz",
  t1_post = "3DT1GADO.nii.gz",
  flair = "3DFLAIR.nii.gz",
  t2 = "T2.nii.gz",
  pd = "DP.nii.gz",
  gold_standard = NULL,
  outdir = ".",
  num_templates = 15,
  verbose = TRUE){

  process_images(t1_pre = t1_pre,
                 t1_post = t1_post,
                 flair = flair,
                 t2 = t2,
                 pd = pd,
                 gold_standard = gold_standard,
                 outdir = outdir,
                 num_templates = num_templates,
                 verbose = verbose)
  # rm(list = "proc_L")
  # L = list(mask = ret_mask_fname,
  #          gold_standard = les_out_fname,
  #          all_imgs =   every_fname,
  #          outdir = outdir)


  have_gold_standard = !is.null(gold_standard)
  L = make_processed_filenames(
    types = c("FLAIR", "T1_Pre", "T1_Post", "T2", "PD"),
    gold_standard = have_gold_standard,
    outdir = outdir)


  df = make_df(fnames = L$fnames,
               mask = L$mask,
               y_fname = L$gold_standard)

  mask = check_nifti(L$mask)

  ind = df$qFLAIR > 0.5

  p = rep(0, length = nrow(df))
  ##################################
  # Predict from the model
  ##################################
  p[ind] = predict(msseg::rf_model,
                   newdata = df[ind,],
                   type = "prob")[, "1"]
  pimg = remake_img(p, img = mask, mask = mask)
  yhat = pimg > (msseg::cutoff)

  yhat = unstrip_image(strip_image = yhat,
                       outdir = outdir)
  return(yhat)
}








