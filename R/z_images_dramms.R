#' @title Standardized to Template Images DRAMMS
#' @description Normalizes to the template
#' @param t1_inv_deffile Name of inverse deformation file
#' @param orig_t1_pre Skull-stripped T1 image to register to Eve
#' @param t1_pre Normalized T1 Pre
#' @param t1_post Normalized T1 Post
#' @param flair Normalized FLAIR
#' @param t2 Normalized T2
#' @param pd Normalized PD
#' @param outdir output directory
#' @param add_stub Stub to add on the filenames
#'
#' @return NULL
#' @export
#' @importFrom fslr checkimg readnii finite_img robust_window writenii
#' @import drammsr
z_images_dramms = function(t1_inv_deffile,
                           orig_t1_pre,
                           t1_pre,
                           t1_post,
                           flair,
                           t2,
                           pd,
                           outdir = ".",
                           add_stub = "_strip_norm"){
  t1_inv_deffile = checkimg(t1_inv_deffile)
  verbose = TRUE

  niis = list(FLAIR = flair,
              T1_Pre = t1_pre,
              T1_Post = t1_post,
              T2 = t2,
              PD = pd)

  # REMOVE NULL
  nulls = sapply(niis, is.null)
  niis = niis[!nulls]

  keep_cols = nii_names = names(niis)

  stopifnot(length(keep_cols) == length(niis))

  imgs = check_nifti(niis)
  stopifnot(all(c("T1_Pre", "FLAIR") %in% nii_names))
  if (verbose) {
    message("Reading in Images")
  }

  #############################################
  # Columns to keep
  #############################################
  ofnames = file.path(outdir, keep_cols)
  ofnames = outer(ofnames,
                  c("_Z_Mean", "_Z_Median"),
                  paste0)
  ofnames = paste0(ofnames,
                   ".nii.gz")


  # if (!all_exists(ofnames))) {
  #############################################
  # Get stubs for whitestripes
  #############################################
  # reg_name = "T1_Pre"
  # reg_file = irow[, reg_name]
  # template_brain = system.file("JHU_MNI_SS_T1_brain.nii.gz",
  #                              package = "msseg")
  #############################################
  # Get stubs for whitestripes
  #############################################
  # outprefix = tempfile()
  # temp = readnii(template_brain)

  reg_file = checkimg(orig_t1_pre)
  # native = check_nifti(orig_t1_pre)


  template_dir = system.file("Normed_Templates",
                             package = "msseg")
  for (icol in keep_cols) {

    out_stub = file.path(template_dir,
                         paste0("Template_", icol, "_DRAMMS"))

    out_name = paste0(out_stub, "_Median")
    med_img = readnii(out_name)

    out_name = paste0(out_stub, "_SD")
    sd_img = readnii(out_name)

    out_name = paste0(out_stub, "_Mean")
    mean_img = readnii(out_name)


    ofname = file.path(outdir, paste0(icol, add_stub))
    nim = check_nifti(imgs[[icol]])

    z_mean = finite_img(
      (nim - mean_img) /
        sd_img )

    z_med = finite_img(
      (nim - med_img) /
        sd_img )

    z_mean_native = dramms_warp(
      z_mean,
      def = t1_inv_deffile,
      template = reg_file,
      retimg = TRUE)
    z_mean_native = robust_window(
      z_mean_native,
      probs = c(0.01, 0.99))
    writenii(z_mean_native,
             filename = paste0(ofname,
                               "_DRAMMS",
                               "_Z_Mean.nii.gz"))

    z_med_native = dramms_warp(
      z_med,
      def = t1_inv_deffile,
      template = reg_file,
      retimg = TRUE)
    z_med_native = robust_window(
      z_med_native,
      probs = c(0.01, 0.99))
    writenii(z_med_native,
             filename = paste0(ofname,
                               "_DRAMMS",
                               "_Z_Median.nii.gz"))

    if (verbose) {
      message(icol)
    }
  }
  return(invisible(NULL))
}