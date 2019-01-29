#' @title MALF tissue Segmentation
#' @description Uses MALF for tissue-class segmentation
#'
#' @param t1 T1-weighted image/filename
#' @param outfile Output filename
#' @param num_templates Number of templates to use for MALF
#' @param interpolator Interpolation used
#' @param typeofTransform Transformation to align the templates to the T1
#' @param func Function used to vote
#' @param keep_regs Keep registrations (for future use),
#' passed to \code{\link{malf}}
#' @param verbose Print diagnostic messages
#'
#' @return Object of class nifti
#' @export
#' @importFrom extrantsr malf
#' @importFrom neurobase check_outfile
malf_tissue_seg = function(
  t1,
  outfile = NULL,
  num_templates = 15,
  interpolator = "NearestNeighbor",
  typeofTransform = "SyN",
  func = "mode",
  keep_regs = TRUE,
  verbose = TRUE){

  root_mass_template_dir = system.file("MASS_Templates",
                                       package = "msseg")
  #######################################
  # Try MALF for Tissues with MASS Templates
  #######################################
  template_files = file.path(root_mass_template_dir,
                             "WithCerebellum_noneck",
                             paste0("Template", 1:num_templates,
                                    ".nii.gz"))
  # template_masks = sub("[.]nii",
  #                      "_str_cbq.nii",
  #                      template_files)
  template_structs = sub("[.]nii",
                         "_Tissue_Classes.nii",
                         template_files)
  template_files = sub("[.]nii",
                       "_strip.nii",
                       template_files)

  t1 = check_nifti(t1)
  outfile = check_outfile(outfile, retimg = TRUE)
  fnames = outfile

  if (all_exists(fnames)) {
    tissue_seg = readnii(fnames)
  } else {
    tissue_seg = malf(
      infile = t1,
      keep_images = FALSE,
      template.images = template_files,
      template.structs = template_structs,
      interpolator = interpolator,
      typeofTransform = typeofTransform,
      keep_regs = keep_regs,
      func = func,
      verbose = verbose)
    if (keep_regs) {
      oimg = tissue_seg$outimg
      writenii(oimg,
               filename = fnames)
    } else {
      writenii(tissue_seg,
               filename = fnames)
    }
  }
  return(tissue_seg)
}
