#' @title MALF tissue Segmentation without Re-registration
#' @description Uses MALF for tissue-class segmentation, but with already-done
#' registrations
#'
#' @param t1 T1-weighted image/filename
#' @param regs List of registrations from malf or malf_registration,
#' each element must have fwdtransforms and interpolator.
#' Length must match \code{num_templates}.
#' @param outfile Output filename
#' @param num_templates Number of templates to use for MALF
#' @param interpolator Interpolation used
#' @param typeofTransform Transformation to align the templates to the T1
#' @param func Function used to vote
#' @param verbose Print diagnostic messages
#' @param ... additional options to pass to \code{reapply_malf}
#'
#' @return Object of class nifti
#' @export
#' @importFrom extrantsr malf
#' @importFrom neurobase check_outfile
#' @importFrom extrantsr reapply_malf
reapply_malf_tissue_seg = function(t1,
                                   regs,
                                   outfile = NULL,
                                   num_templates = 15,
                                   interpolator = "NearestNeighbor",
                                   typeofTransform = "SyN",
                                   func = "mode",
                                   verbose = TRUE,
                                   ...){

  root_mass_template_dir = system.file("MASS_Templates",
                                       package = "msseg")
  #######################################
  # Try MALF for Tissues with MASS Templates
  #######################################
  template_files = file.path(root_mass_template_dir,
                             "WithCerebellum_noneck",
                             paste0("Template", 1:num_templates,
                                    ".nii.gz"))
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
    tissue_seg = reapply_malf(
      infile = t1,
      regs = regs,
      template.structs = template_structs,
      keep_images = FALSE,
      outfile = fnames,
      retimg = TRUE,
      func = func,
      interpolator = interpolator,
      verbose = verbose,
      ...
    )
  }
  return(tissue_seg)
}
