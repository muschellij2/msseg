#' @title Process Images for MSSeg
#' @description Processing script for MSSeg Challenge
#' @param t1_pre T1-weighted pre-gad image/filename
#' @param t1_post T1-weighted post-gad image/filename
#' @param flair FLAIR image/filename
#' @param t2 T2-weighted image/filename
#' @param pd Proton Density pre-gad image/filename
#' @param gold_standard Gold Standar Lesion
#' @param outdir Output directory
#' @param num_templates Number of templates used in MASS
#' @param verbose print diagnostic messages
#' @param force Force tissue segmentations
#' @return List of output filenames
#' @export
#' @importFrom plyr llply alply
#' @importFrom extrantsr oro2ants ants2oro double_remove_neck bias_correct
#' @importFrom extrantsr oMath filler within_visit_registration reg_flip otropos
#' @importFrom extrantsr corr_img diff_self create_moment
#' @importFrom ANTsR getMask "%>%"
#' @importFrom fslr dropEmptyImageDimensions fslbet quantile_img
#' @importFrom fslr fslsmooth nii.stub remake_img datatyper
process_images = function(t1_pre,
                          t1_post,
                          flair,
                          t2,
                          pd,
                          gold_standard = NULL,
                          outdir = ".",
                          num_templates = 15,
                          verbose = TRUE,
                          force = TRUE){
  verbose = TRUE

  niis = list(FLAIR = flair,
              T1_Pre = t1_pre,
              T1_Post = t1_post,
              T2 = t2,
              PD = pd)

  # REMOVE NULL
  nulls = sapply(niis, is.null)
  niis = niis[!nulls]

  nii_names = names(niis)

  stopifnot(all(c("T1_Pre", "FLAIR") %in% nii_names))
  if (verbose) {
    message("Reading in Images")
  }
  imgs = llply(niis, check_nifti,
               .progress = "text")

  have_gold_standard = !is.null(gold_standard)
  if (have_gold_standard) {
    if (verbose > 0) {
      message("Reading in Gold Standard")
    }
    les = check_nifti(gold_standard)
  }
  rm(list = "niis")

  if (verbose > 0) {
    message("N4 Bias Field Correction")
  }
  #################################
  # Bias correct the images using N4
  #################################
  fnames = file.path(outdir,
                     paste0(nii_names,
                            "_N4.nii.gz"))
  names(fnames) = nii_names

  if (all_exists(fnames)) {
    n4 = llply(fnames, readnii,
               .progress = "text")
  } else {
    n4 = llply(imgs,
               bias_correct,
               correction = "N4",
               verbose = verbose,
               .progress = "text")

    #################################
    # Swap image for remove_neck.
    # Can use swapdim in remove_neck, but
    # registration seems to work better if
    # images are in std RPI
    #################################
    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, n4, fnames)
  }
  names(n4) = nii_names

  rm(list = "imgs")

  if (verbose > 0) {
    message("Removing necks")
  }
  ####################################################
  # Making fname stubs
  ####################################################
  fname_stubs = file.path(outdir,
                          names(n4))
  names(fname_stubs) = names(n4)

  ####################################################
  # Removing Neck
  ####################################################
  fnames = paste0(fname_stubs,
                  "_noneck.nii.gz")
  names(fnames) = names(n4)

  if (all_exists(fnames)) {
    noneck = llply(fnames,
                   readnii,
                   .progress = "text")
  } else {


    mni.template.file = system.file("MNI152_T1_1mm_brain.nii.gz",
                                    package = "msseg")
    mni.template.mask = system.file("MNI152_T1_1mm_brain_mask.nii.gz",
                                    package = "msseg")


    noneck = llply(n4,
                   double_remove_neck,
                   template.file = mni.template.file,
                   template.mask = mni.template.mask,
                   .progress = "text")
    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, noneck, fnames)
  }

  if (verbose > 0) {
    message("Dropping Empty Dimensions")
  }
  ####################################################
  # Dropping empty dimensions
  ####################################################
  fnames = paste0(fname_stubs,
                  "_reduced.nii.gz")
  names(fnames) = names(n4)

  if (all_exists(fnames)) {
    rm_neck = llply(fnames,
                    readnii,
                    .progress = "text")
  } else {
    drops = llply(noneck,
                  function(nn){
                    x = oro2ants(nn)
                    gm = getMask(x, cleanup = 0)
                    x = ants2oro(gm)
                    dd = dropEmptyImageDimensions(
                      x, keep_ind = TRUE,
                      other.imgs = nn)
                    dd$outimg = dd$other.imgs
                    dd$other.imgs = NULL
                    dd
                  },
                  .progress = "text")
    rm_neck = llply(drops,
                    function(x){
                      x$outimg
                    })

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, rm_neck, fnames)
  }

  #################################
  # Removing dims if gold standard exists
  #################################
  if (have_gold_standard) {
    les_fname = file.path(outdir,
                          "GoldStandard_reduced.nii.gz")
    if (all_exists(les_fname)) {
      les_rm_neck = readnii(les_fname)
    } else {
      nn = noneck$FLAIR
      x = oro2ants(nn)
      gm = getMask(x, cleanup = 0)
      x = ants2oro(gm)

      dd = dropEmptyImageDimensions(
        x, other.imgs = les)
      les_rm_neck = dd$other.imgs
      writenii(les_rm_neck,
               filename = les_fname)
    }
  }
  rm(list = "noneck"); gc(); gc();

  if (verbose > 0) {
    message("Getting Quick mask")
  }
  #######################################
  # Just see how a quick mask can be created
  # with getMask and FLAIR
  #######################################
  fnames = file.path(outdir,
                     "FLAIR_reduced_quickmask.nii.gz")
  if (all_exists(fnames)){
    quick_flair_mask = readnii(fnames)
  } else {
    quick_flair_mask = oro2ants(rm_neck$FLAIR)
    quick_flair_mask = getMask(quick_flair_mask)
    quick_flair_mask = ants2oro(quick_flair_mask)
    writenii(quick_flair_mask,
             filename = fnames)
  }

  #################################
  # Filling
  #################################
  fnames = file.path(outdir,
                     "FLAIR_reduced_quickmask_filled.nii.gz")

  if (all_exists(fnames)) {
    filled_quick_flair_mask = readnii(fnames)
  } else {
    filled_quick_flair_mask = quick_flair_mask %>%
      oMath("MD", 5) %>% oMath("ME", 5)
    writenii(filled_quick_flair_mask,
             filename = fnames)
  }


  #################################
  # Running BET
  #################################
  fnames = file.path(outdir,
                     "FLAIR_BET.nii.gz")
  if (all_exists(fnames)) {
    # flair_bet = readnii(fnames)
    flair_bet = fnames
  } else {
    flair_bet = fslbet(rm_neck$FLAIR)
    writenii(flair_bet,
             filename = fnames)
  }
  rm(list = "flair_bet"); gc(); gc();

  #######################################
  # Need MASS Templates
  #######################################
  # root_mass_template_dir = system.file("MASS_Templates",
  #                                      package = "msseg")
  # mass_template_dir = file.path(root_mass_template_dir,
  #                               "WithCerebellum")
  # template_files = file.path(
  #   mass_template_dir,
  #   paste0("Template", 1:num_templates,
  #          ".nii.gz"))
  # template_structs = sub("[.]nii",
  #                        "_str_cbq.nii",
  #                        template_files)
  #
  # #######################################
  # # Try MALF with MASS Templates
  # #######################################
  # fnames = paste0(fname_stubs,
  #                 "_MALF_Brain_Mask.nii.gz")
  # names(fnames) = names(n4)
  #
  # mask_fnames = fnames
  # if (!all_exists(fnames)) {
  #
  #   keep = !file.exists(fnames)
  #   run_rm_neck = rm_neck[keep]
  #   run_fnames = fnames[keep]
  #   brain_masks = vector(mode = "list",
  #                        length = length(run_fnames))
  #   names(brain_masks) = names(run_fnames)
  #
  #   for (i in seq_along(run_fnames)) {
  #     # brain_masks = llply(run_rm_neck,
  #     brain_masks[[i]] = malf(
  #       run_rm_neck[[i]],
  #       keep_images = FALSE,
  #       template.images =
  #         template_files,
  #       template.structs =
  #         template_structs,
  #       interpolator =
  #         "NearestNeighbor",
  #       typeofTransform = "SyN")
  #     writenii(brain_masks[[i]],
  #              filename = run_fnames[i])
  #   }
  # }

  if (verbose > 0) {
    message("Running BET to get Brain mask for each image")
  }
  #######################################
  # Try BET
  #######################################
  fnames = paste0(fname_stubs,
                  "_MALF_Brain_Mask.nii.gz")
  names(fnames) = names(n4)

  mask_fnames = fnames
  if (!all_exists(fnames)) {

    keep = !file.exists(fnames)
    run_rm_neck = rm_neck[keep]
    run_fnames = fnames[keep]
    brain_masks = vector(mode = "list",
                         length = length(run_fnames))
    names(brain_masks) = names(run_fnames)

    for (i in seq_along(run_fnames)){
      # brain_masks = llply(run_rm_neck,
      bet_mask = fslbet(run_rm_neck[[i]]) > 0
      bet_mask = filler(bet_mask, fill_size = 2)
      brain_masks[[i]] = bet_mask
      writenii(brain_masks[[i]],
               filename = run_fnames[i])
      rm(list = "bet_mask"); gc(); gc();
    }
    rm(list = "run_rm_neck"); gc(); gc();
  }

  #######################################
  # Read in Brain Masks
  #######################################
  if (all_exists(mask_fnames)){
    brain_masks = llply(mask_fnames,
                        readnii,
                        .progress = "text")
  } else {
    stop("Need masks!")
  }

  if (verbose > 0) {
    message("Running N4 on Brains only")
  }
  ###################################################
  # Bias correct the images using N4, only on the brain
  ###################################################
  fnames = paste0(fname_stubs,
                  "_Brain_N4.nii.gz")
  names(fnames) = names(n4)

  if (all_exists(fnames)){
    n4_brain = llply(fnames,
                     readnii,
                     .progress = "text")
  } else {
    n4_brain = mapply(function(img,
                               mask) {
      bias_correct(
        file = img,
        mask = mask,
        correction = "N4"
      )
    }, rm_neck, brain_masks,
    SIMPLIFY = FALSE)

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, n4_brain, fnames)
  }


  if (verbose > 0) {
    message("Running Rigid Registration to FLAIR")
  }
  ################################
  # Register to FLAIR Image
  ################################
  to_img = "FLAIR"
  from_imgs = c("T1_Pre", "T1_Post",
                "T2", "PD")
  from_imgs = intersect(from_imgs, nii_names )
  fixed = n4_brain[[to_img]]
  mask_full = brain_masks[[to_img]]
  img_names = c(to_img, from_imgs)

  every_fname = NULL
  fnames = file.path(outdir,
                     paste0(img_names,
                            "_strip.nii.gz"))
  names(fnames) = img_names
  every_fname = c(every_fname, fnames)
  mask_fname = file.path(outdir,
                         "FLAIR_MALF_Brain_Mask_strip.nii.gz")
  ret_mask_fname = mask_fname
  if (all_exists(c(
    fnames, mask_fname))) {
    masked_reg_imgs = llply(fnames,
                            readnii,
                            .progress = "text")
    mask = readnii(mask_fname)
  } else {
    reg = within_visit_registration(
      fixed = fixed,
      moving = n4_brain[from_imgs],
      typeofTransform = "Rigid",
      interpolator =
        "LanczosWindowedSinc",
      outfiles = NULL)

    reg_imgs = lapply(reg, `[[`,
                      "outfile")
    reg_imgs = c(list(fixed), reg_imgs)
    names(reg_imgs) = img_names

    masked_reg_imgs = llply(reg_imgs,
                            mask_img, mask = mask_full,
                            .progress = "text")

    #############################
    # Should we drop again here
    # for brain_masks?
    #############################
    dd_all = dropEmptyImageDimensions(
      mask_full,
      other.imgs =
        masked_reg_imgs)

    masked_reg_imgs = dd_all$other.imgs

    mask = dd_all$outimg
    writenii(mask, filename = mask_fname)

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, masked_reg_imgs, fnames)
    rm(list = "dd_all"); gc();
  }
  rm(list = "n4_brain"); gc(); gc();

  #############################################
  # Dropping dimensions from lesion
  #############################################
  if (have_gold_standard) {
    les_fname = file.path(outdir,
                          "GoldStandard_strip.nii.gz")
    les_out_fname = les_fname
    if (all_exists(les_fname)) {
      # les_strip = readnii(les_fname)
      les_strip = les_fname
    } else {
      dd = dropEmptyImageDimensions(
        mask_full,
        other.imgs = les_rm_neck)

      les_strip = dd$other.imgs
      writenii(les_strip,
               filename = les_fname)
      rm(list = "les_strip"); gc();
      les_strip = les_fname
    }
  }
  rm(list = "rm_neck"); gc(); gc();
  if (verbose > 0) {
    message("Running Tissue Segmentation MALF on T1 Pre")
  }
  #######################################
  # Try MALF for Tissues with MASS Templates
  #######################################
  fnames = file.path(outdir,
                     "T1_Pre_MALF_Tissue_Classes.nii.gz")
  every_fname = c(every_fname, fnames)

  tissue_seg_reg = malf_tissue_seg(
    t1 = masked_reg_imgs$T1_Pre,
    outfile = fnames,
    num_templates = num_templates,
    interpolator = "NearestNeighbor",
    typeofTransform = "SyN",
    func = "mode",
    keep_regs = TRUE,
    verbose = TRUE)
  # tissue_seg = tissue_seg_reg$outimg

  if (verbose > 0) {
    message("Running Tissue Segmentation MALF on T1 Pre with Gauss")
  }
  #######################################
  # Try MALF for Tissues with MASS Templates
  #######################################
  fnames = file.path(outdir,
                     "T1_Pre_MALF_Tissue_Classes_Gauss.nii.gz")
  every_fname = c(every_fname, fnames)

  if (inherits(tissue_seg_reg, "nifti") || force) {
    regs = list(fwdtransforms = tempfile())
  } else {
    regs = tissue_seg_reg$regs
  }
  rm(list = "tissue_seg_reg"); gc(); gc();
  trans = regs$fwdtransforms
  #########################################
  # Saves computation by not having to redo registrations
  #########################################
  if (all_exists(trans) && !force) {
    tissue_seg_gauss = reapply_malf_tissue_seg(
      t1 = masked_reg_imgs$T1_Pre,
      regs = regs,
      outfile = fnames,
      num_templates = num_templates,
      interpolator = "MultiLabel",
      func = "mode",
      verbose = TRUE)
  } else {
    tissue_seg_gauss = malf_tissue_seg(
      t1 = masked_reg_imgs$T1_Pre,
      outfile = fnames,
      num_templates = num_templates,
      interpolator = "MultiLabel",
      typeofTransform = "SyN",
      func = "mode",
      keep_regs = FALSE,
      verbose = TRUE)
  }
  rm(list = "tissue_seg_gauss"); gc(); gc();
  if (verbose > 0) {
    message("Running Atropos Tissue Segmentation on T1 Pre")
  }
  #######################################
  # Try Atropos
  #######################################
  fnames = file.path(outdir,
                     "T1_Pre_ANTs_Tissue_Classes.nii.gz")
  every_fname = c(every_fname, fnames)

  tissue_classes = 1:3
  tissue_fnames = file.path(outdir,
                            paste0("T1_Pre_ANTs_Tissue_Prob_",
                                   tissue_classes, ".nii.gz") )
  names(tissue_fnames) = paste0("prob_",
                                tissue_classes)
  every_fname = c(every_fname, tissue_fnames)


  if (all_exists(c(fnames, tissue_fnames))) {
    # tissue_ants = readnii(fnames)
    tissue_ants = fnames
    tissue_probs = llply(tissue_fnames,
                         readnii,
                         .progress = "text")
  } else {
    rr = robust_window(masked_reg_imgs$T1_Pre)
    tissue_ants = otropos(rr,
                          x = mask,
                          v = 1)
    tissue_probs =
      tissue_ants$probabilityimages
    tissue_ants = tissue_ants$segmentation
    writenii(tissue_ants,
             filename = fnames)

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, tissue_probs, tissue_fnames)
    rm(list = c("tissue_ants", "rr")); gc()
    tissue_ants = fnames
  }

  if (verbose > 0) {
    message("Running Atropos Tissue Segmentation on FLAIR")
  }
  ##################################
  # Use FLAIR Dummy
  ##################################
  fnames = file.path(outdir,
                     "FLAIR_ANTs_Tissue_Classes.nii.gz")
  every_fname = c(every_fname, fnames)

  flair_tissue_classes = 1:4
  flair_tissue_fnames = file.path(
    outdir,
    paste0("FLAIR_ANTs_Tissue_Prob_",
           flair_tissue_classes, ".nii.gz") )
  names(flair_tissue_fnames) = paste0(
    "flair_prob_",
    flair_tissue_classes)
  every_fname = c(every_fname, flair_tissue_fnames)



  if (all_exists(c(fnames,
                   flair_tissue_fnames))) {
    # flair_tissue_ants = readnii(fnames)
    flair_tissue_ants = fnames
    flair_tissue_probs = llply(
      flair_tissue_fnames,
      readnii,
      .progress = "text")
  } else {
    rr = robust_window(
      masked_reg_imgs$FLAIR)
    flair_tissue_ants = otropos(rr,
                                x = mask,
                                i = paste0("KMeans[",
                                           length(flair_tissue_classes),
                                           "]"),
                                v = 1)
    flair_tissue_probs =
      flair_tissue_ants$probabilityimages
    flair_tissue_ants = flair_tissue_ants$segmentation
    writenii(flair_tissue_ants,
             filename = fnames)

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, flair_tissue_probs,
    flair_tissue_fnames)
    rm(list = c("flair_tissue_ants", "rr")); gc()
    flair_tissue_ants = fnames
  }

  if (verbose > 0) {
    message("Running Atropos Tissue Segmentation on T1 Post")
  }
  ##################################
  # Use T1 Post
  ##################################
  fnames = file.path(outdir,
                     "T1_Post_ANTs_Tissue_Classes.nii.gz")
  every_fname = c(every_fname, fnames)

  if (all_exists(fnames)) {
    # post_tissue_ants = readnii(fnames)
  } else {
    rr = robust_window(
      masked_reg_imgs$T1_Post)
    post_tissue_ants = otropos(rr,
                               i = "KMeans[4]",
                               x = mask,
                               v = 1)
    post_tissue_ants =
      post_tissue_ants$segmentation
    writenii(post_tissue_ants,
             filename = fnames)
    rm(list = c("post_tissue_ants", "rr")); gc()
  }

  if (verbose > 0) {
    message("Quantiling Image")
  }
  ################################
  # Quantile Images
  ################################
  fnames = file.path(outdir,
                     paste0(img_names,
                            "_strip_quantile.nii.gz"))
  names(fnames) = img_names
  every_fname = c(every_fname, fnames)


  if (all_exists(c(fnames)))  {
    # qimgs = llply(fnames,
    #               readnii,
    #               .progress = "text")
  } else {
    qimgs = llply(masked_reg_imgs,
                  quantile_img,
                  mask = mask,
                  .progress = "text")
    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, qimgs, fnames)
  }
  qimgs = fnames

  ################################
  # Delete CSF from Mask
  ################################
  fnames = file.path(outdir,
                     paste0("Over_FLAIR_15Quantile.nii.gz"))
  flair = masked_reg_imgs$FLAIR
  f90_fname = file.path(outdir,
                        "FLAIR_90.nii.gz")
  every_fname = c(every_fname, fnames, f90_fname)


  if (all_exists(c(fnames, f90_fname)))  {
    mask2 = readnii(fnames)
    thresh = readnii(f90_fname)
  } else {
    vals = flair[mask == 1]
    q15 = quantile(vals,
                   probs = c(0.15, 0.9))
    mask2 = mask &
      flair > q15[1]
    mask2 = datatyper(mask2)
    mask2 = mask2 %>%
      oMath("GetLargestComponent")

    writenii(mask2,
             filename = fnames
    )

    thresh = flair > q15[2]
    thresh = datatyper(thresh)
    thresh = mask_img(thresh,
                      mask)
    writenii(thresh,
             filename = f90_fname
    )
    # don't delete mas2
    rm(list = c("vals", "thresh")); gc(); gc();
  }

  if (verbose > 0) {
    message("Normalizing Images")
  }
  ################################
  # Normalize Using Trimmed Brain Mask
  ################################
  fnames = file.path(outdir,
                     paste0(img_names,
                            "_strip_norm.nii.gz"))
  names(fnames) = img_names
  norm_fnames = fnames
  every_fname = c(every_fname, fnames)


  if (all_exists(fnames)) {
    norm_imgs = llply(fnames,
                      readnii,
                      .progress = "text")
  } else {
    norm_imgs = llply(
      masked_reg_imgs,
      trimmed_z,
      mask = mask2,
      trim = 0.2,
      .progress = "text")

    norm_imgs = llply(
      norm_imgs,
      mask_img,
      mask = mask,
      .progress = "text")

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, norm_imgs, fnames)

  }

  all_imgs = norm_imgs
  names(all_imgs) = names(norm_imgs)

  all_fnames = c(
    norm_fnames)
  all_fname_stubs = nii.stub(all_fnames)
  all_fnames = paste0(
    all_fname_stubs,
    "_toEve.nii.gz")
  names(all_fnames) = names(all_imgs)
  names(all_fname_stubs) = names(all_imgs)

  # ################################
  # # Register to Eve Template
  # Nonlinearly with T1
  # ################################
  # template_brain = file.path(
  #   eve_template_dir,
  #   "JHU_MNI_SS_T1_brain.nii.gz")
  if (verbose > 0) {
    message("Registering to Eve")
  }
  template_brain = system.file("JHU_MNI_SS_T1_brain.nii.gz",
                               package = "msseg")

  eve_fnames = fnames = all_fnames

  reg_name = "T1_Pre"
  other.names = names(all_fnames)
  other.names = setdiff(other.names,
                        reg_name)
  t1 = masked_reg_imgs[[reg_name]]

  if (!all_exists(fnames)) {
    # t1 = masked_reg_imgs[[reg_name]]
    t1_reg = registration(
      filename = t1,
      template.file = template_brain,
      typeofTransform = "SyN",
      interpolator = "LanczosWindowedSinc",
      # outfile = fnames[reg_name],
      outfile = tempfile(fileext = ".nii.gz"),
      other.files =
        all_imgs[other.names],
      other.outfiles =
        fnames[other.names],
      remove.warp = FALSE,
      outprefix = tempfile()
    )
    rm(t1_reg); gc(); gc()
  }

  if (verbose > 0) {
    message("Z-scoring from Eve")
  }
  fnames = nii.stub(norm_fnames)
  fnames = outer(fnames,
                 c("_Z_Mean", "_Z_Median"),
                 paste0)
  fnames = paste0(fnames,
                  ".nii.gz")
  every_fname = c(every_fname, fnames)

  # ################################
  # Inverting the T1 from DRAMMS
  # ################################
  if (!all_exists(fnames)) {
    z_images(orig_t1_pre = t1,
             t1_pre = eve_fnames["T1_Pre"],
             t1_post = eve_fnames["T1_Post"],
             flair = eve_fnames["FLAIR"],
             t2 = eve_fnames["T2"],
             pd = eve_fnames["PD"],
             outdir = outdir,
             add_stub = "_strip_norm")
  }


  # ################################
  # # Register to Eve Template
  # DRAMMS with T1
  # ################################
  if (verbose > 0) {
    message("DRAMMS to Eve")
  }
  dramms_fnames =
    fnames =
    gsub("toEve","toEve_DRAMMS", all_fnames)
  reg_name = "T1_Pre"
  other.names = names(all_fnames)
  other.names = setdiff(other.names,
                        reg_name)
  t1_outfile = fnames[reg_name]
  t1_deffile = paste0(
    nii.stub(t1_outfile),
    "_def.nii.gz")

  if (!all_exists(fnames)) {

    d = dramms(
      source = t1,
      target = template_brain,
      outfile = t1_outfile,
      outdef = t1_deffile
    )
    mapply(function(img, outfile){
      dramms_warp(source = img,
                  def = t1_deffile,
                  outfile = outfile,
                  template = template_brain)
    },  all_imgs[other.names],
    fnames[other.names])
  }

  # ################################
  # Register to T1
  # DRAMMS with Eve Template
  # ################################
  if (verbose > 0) {
    message("DRAMMS Eve to T1")
  }
  t1_inv_deffile = paste0(
    nii.stub(t1_outfile),
    "_inv_def.nii.gz")
  template_to_pre = file.path(
    dirname(t1_outfile),
    "Eve_to_T1.nii.gz"
  )
  if (!all_exists(c(t1_inv_deffile,
                    template_to_pre))) {
    d = dramms(
      source = template_brain,
      target = t1,
      outfile = template_to_pre,
      outdef = t1_inv_deffile
    )
  }

  if (verbose > 0) {
    message("Inverting DRAMMS Z-score to Eve")
  }
  # ################################
  # Inverting the T1 from DRAMMS
  # ################################
  t1_inv_deffile = paste0(
    nii.stub(t1_outfile),
    "_inv_def.nii.gz")

  fnames = nii.stub(norm_fnames)
  fnames = outer(fnames,
                 c("_DRAMMS_Z_Mean", "_DRAMMS_Z_Median"),
                 paste0)
  fnames = paste0(fnames,
                  ".nii.gz")
  every_fname = c(every_fname, fnames)

  if (!all_exists(fnames)) {
    z_images_dramms(t1_inv_deffile = t1_inv_deffile,
                    orig_t1_pre = t1,
                    t1_pre = dramms_fnames["T1_Pre"],
                    t1_post = dramms_fnames["T1_Post"],
                    flair = dramms_fnames["FLAIR"],
                    t2 = dramms_fnames["T2"],
                    pd = dramms_fnames["PD"],
                    outdir = outdir,
                    add_stub = "_strip_norm")
  }




  if (verbose > 0) {
    message("RAVENS")
  }
  reg_name = "T1_Pre"
  template_mask = readnii(template_brain) > 0

  template_ravens_pre = file.path(
    outdir,
    "Eve_to_T1_RAVENS"
  )
  fnames = paste0(
    template_ravens_pre,
    ".nii.gz")
  if (!all_exists(fnames)) {
    t1 = masked_reg_imgs[[reg_name]]
    d = dramms_with_ravens(
      source = template_brain,
      target = t1,
      outfile = NULL,
      outdef = tempfile(fileext =
                          ".nii.gz"),
      label = template_mask,
      labs = 1,
      ravens_prefix = template_ravens_pre,
      retimg = FALSE
    )
  }

  t1_inv_deffile = paste0(
    nii.stub(t1_outfile),
    "_inv_def.nii.gz")
  template_to_pre = file.path(
    dirname(t1_outfile),
    "Eve_to_T1.nii.gz"
  )
  if (!all_exists(c(t1_inv_deffile,
                    template_to_pre))) {
    reg_name = "T1_Pre"
    t1 = masked_reg_imgs[[reg_name]]
    d = dramms(
      source = template_brain,
      target = t1,
      outfile = template_to_pre,
      outdef = t1_inv_deffile
    )
  }

  # eve_norm_imgs = llply(
  #     fnames,
  #     readnii,
  #     .progress = "text")

  if (verbose > 0) {
    message("Template Flipped Difference")
  }
  # ################################
  # Flipped Difference
  # ################################
  fnames = paste0(all_fname_stubs,
                  "_Flipped_Difference.nii.gz")
  names(fnames) = names(all_fname_stubs)
  every_fname = c(every_fname, fnames)

  if (all_exists(fnames)) {
    # flip_imgs = llply(fnames,
    #                   readnii,
    #                   .progress = "text")
    flip_imgs = fnames
  } else {

    reg_name = "T1_Pre"
    other.names = names(fnames)
    other.names = setdiff(other.names,
                          reg_name)
    t1 = masked_reg_imgs[[reg_name]]
    t1_freg = reg_flip(
      t1 = t1,
      register = TRUE,
      native = TRUE,
      template.file = template_brain,
      typeofTransform = "Affine",
      interpolator = "LanczosWindowedSinc",
      # t1.outfile = fnames[reg_name],
      t1.outfile = tempfile(fileext = ".nii.gz"),
      other.files =
        all_imgs,
      other.outfiles =
        fnames,
      a = "-x",
      b = "y",
      c = "z"
    )
    rm(t1_freg); gc(); gc()
    flipped = llply(fnames,
                    readnii,
                    .progress = "text")

    ##########################
    # Take difference
    ##########################
    flip_imgs = mapply(function(res,
                                fres){
      res - fres
    }, all_imgs, flipped,
    SIMPLIFY = FALSE)
    rm(list = "flipped"); gc(); gc();

    ##########################
    # Mask flipped difference image
    ##########################
    flip_imgs = lapply(flip_imgs,
                       mask_img,
                       mask = mask)

    flip_imgs = flip_imgs[names(fnames)]
    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, flip_imgs, fnames)
    rm(list = "flip_imgs"); gc();
  }

  rm(list = "masked_reg_imgs"); gc(); gc();

  ##############################
  # Smooth Image with 3mm
  # Gaussian and write out
  ##############################
  if (verbose > 0) {
    message("Running Smoothers")
  }
  sigma = 3
  sigmas = c(3, 10, 20)

  # tmp_list = all_imgs
  # tmp_list = lapply(tmp_list, function(x){
  #   NULL
  # })
  sigma_list = vector(
    mode = "list",
    length = length(sigmas))
  isigma = 1

  for (isigma in seq_along(sigmas)) {
    sigma = sigmas[isigma]
    fnames = paste0(all_fname_stubs,
                    "_smooth", sigma, ".nii.gz")
    names(fnames) = names(all_fname_stubs)
    every_fname = c(every_fname, fnames)

    if (all_exists(fnames)){
      # sig_imgs = llply(fnames,
      #                  readnii,
      #                  .progress = "text")
    } else {
      sig_imgs = llply(all_imgs,
                       fslsmooth,
                       sigma = sigma,
                       mask = mask,
                       .progress = "text")

      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, sig_imgs, fnames)
      rm(list = "sig_imgs"); gc()
    }
    sigma_list[[isigma]] = fnames
    # sig_imgs
    print(sigma)
  }


  ##############################
  # Perona Malik Smoothers
  ##############################
  if (verbose > 0) {
    message("Running Perona Malik Smoother")
  }
  fnames = paste0(all_fname_stubs,
                  "_PeronaMalik.nii.gz")
  names(fnames) = names(all_fname_stubs)
  every_fname = c(every_fname, fnames)

  if (all_exists(fnames)) {
    # pm_imgs = llply(fnames,
    #                 readnii,
    #                 .progress = "text")
  } else {

    ##########################
    # Smooth the images
    ##########################
    pm_imgs = llply(all_imgs,
                    oMath,
                    operation = "PeronaMalik",
                    30, 2,
                    .progress = "text")

    pm_imgs = pm_imgs[names(fnames)]
    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, pm_imgs, fnames)
    rm(list = "pm_imgs"); gc();
  }


  ##############################
  # Creating Moment Images
  ##############################
  if (verbose > 0) {
    message("Creating Moment Images")
  }
  ifname = 1
  all_mom_imgs = vector(
    mode = "list",
    length = length(all_fname_stubs))
  names(all_mom_imgs) = names(all_fname_stubs)

  for (ifname in seq_along(all_fname_stubs)) {

    fname_stub = all_fname_stubs[ifname]
    stubs = c("mn", "sd", "skew",
              "kurt", "grad", "z")
    fnames = paste0(fname_stub, "_",
                    stubs, ".nii.gz")
    names(fnames) = stubs
    every_fname = c(every_fname, fnames)

    if (all_exists(fnames)) {
      # keep this
      mom_imgs = llply(fnames,
                       readnii,
                       .progress = "text")
    } else {
      img = all_imgs[[ifname]]
      mom_imgs = create_moment(img,
                               mask = mask,
                               retimg = TRUE)
      mom_imgs = mom_imgs[stubs]
      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, mom_imgs, fnames)
    }
    all_mom_imgs[[ifname]] = mom_imgs

    print(names(all_mom_imgs)[ifname])
  }

  ################################
  # Quantile Images
  ################################
  mean_imgs = lapply(all_mom_imgs,
                     function(x){
                       x$mn
                     })
  rm(list = "all_mom_imgs"); gc()
  fnames = paste0(nii.stub(norm_fnames),
                  "_mn_quantile.nii.gz")
  names(fnames) = img_names
  every_fname = c(every_fname, fnames)

  if (all_exists(c(fnames))) {
    # mn_qimgs = llply(fnames,
    #                  readnii,
    #                  .progress = "text")
  } else {
    mn_qimgs = llply(mean_imgs,
                     quantile_img,
                     mask = mask,
                     .progress = "text")
    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, mn_qimgs, fnames)
    rm(list = "mn_qimgs"); gc()
  }
  mn_qimgs = fnames
  rm(list = "mean_imgs"); gc();

  ###########################################
  # Creating Moment Images of
  # Probability Class
  ###########################################
  if (verbose > 0) {
    message("Creating Moment Images of Probability Images")
  }
  all_prob_mom_imgs = vector(
    mode = "list",
    length = length(tissue_classes))
  names(all_prob_mom_imgs) = paste0("prob_",
                                    tissue_classes)


  ifname = 1

  for (ifname in seq_along(
    all_prob_mom_imgs)) {

    fname_stub = nii.stub(
      tissue_fnames[ifname])
    stubs = c("mn", "sd", "skew",
              "kurt", "grad", "z")
    fnames = paste0(fname_stub, "_",
                    stubs, ".nii.gz")
    names(fnames) = stubs
    every_fname = c(every_fname, fnames)

    if (all_exists(fnames)) {
      # prob_mom_imgs = llply(fnames,
      #                       readnii,
      #                       .progress = "text")
    } else {
      img = tissue_probs[[ifname]]
      prob_mom_imgs = create_moment(img,
                                    mask = mask,
                                    retimg = TRUE)
      rm(list = "img"); gc();
      prob_mom_imgs = prob_mom_imgs[stubs]
      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, prob_mom_imgs, fnames)
      rm(list = "prob_mom_imgs"); gc()
    }
    all_prob_mom_imgs[[ifname]] =
      fnames

    print(names(all_prob_mom_imgs)[ifname])
  }
  rm(list = "tissue_probs"); gc(); gc();



  ###########################################
  # Creating Moment Images of
  # Probability Class
  ###########################################
  if (verbose > 0) {
    message("Creating Moment Images of FLAIR Probability Images")
  }
  all_flair_prob_mom_imgs = vector(
    mode = "list",
    length = length(flair_tissue_classes))
  names(all_flair_prob_mom_imgs) = paste0(
    "flair_prob_",
    flair_tissue_classes)


  ifname = 1

  for (ifname in seq_along(
    all_flair_prob_mom_imgs)) {

    fname_stub = nii.stub(
      flair_tissue_fnames[ifname])
    stubs = c("mn", "sd", "skew",
              "kurt", "grad", "z")
    fnames = paste0(fname_stub, "_",
                    stubs, ".nii.gz")
    names(fnames) = stubs
    every_fname = c(every_fname, fnames)

    if (all_exists(fnames)) {
      # prob_mom_imgs = llply(fnames,
      #                       readnii,
      #                       .progress = "text")
    } else {
      img = flair_tissue_probs[[ifname]]
      prob_mom_imgs = create_moment(img,
                                    mask = mask,
                                    retimg = TRUE)
      rm(list = "img"); gc(); gc();
      prob_mom_imgs = prob_mom_imgs[stubs]
      mapply(function(img, fname){
        writenii(img, filename = fname)
      }, prob_mom_imgs, fnames)
      rm(list = "prob_mom_imgs"); gc()
    }
    all_flair_prob_mom_imgs[[ifname]] =
      fnames

    print(names(all_flair_prob_mom_imgs)[ifname])
  }
  rm(list = "flair_tissue_probs"); gc(); gc();



  ##########################################
  # Correlation Images
  ##########################################
  if (verbose > 0) {
    message("Correlation Images")
  }
  all_names = names(norm_imgs)
  eg = expand.grid(img1 = "FLAIR",
                   img2 = setdiff(all_names, "FLAIR"))
  eg = as.matrix(eg)
  colnames(eg) = c("img1", "img2")

  df = as.data.frame(eg,
                     stringsAsFactors = FALSE)
  df$name = paste0(df$img1, "_", df$img2,
                   "_Correlation.nii.gz")
  fnames = file.path(outdir,
                     df$name)
  every_fname = c(every_fname, fnames)

  if (!all_exists(fnames)) {
    res = alply(eg, 1,
                function(iname){
                  img1 = norm_imgs[[iname[1]]]
                  img2 = norm_imgs[[iname[2]]]
                  myres = extrantsr::corr_img(img1, img2,
                                      mask = mask,
                                      radius = rep(1,3),
                                      method = "pearson")
                  gc(); gc();
                  return(myres)
                }, .progress = "text")

    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, res, fnames)
    rm(list = "res"); gc(); gc()
  }

  if (verbose > 0) {
    message("Self Flipping")
  }
  ###########################################
  # Self Flipping
  ###########################################
  fnames = file.path(outdir,
                     df$name)
  fnames = paste0(nii.stub(norm_fnames),
                  "_self_flip.nii.gz")
  every_fname = c(every_fname, fnames)

  if (!all_exists(fnames)) {

    res = llply(norm_imgs,
                function(x) {
                  adder = 5
                  x = mask_img(x + adder,
                               mask)
                  r = extrantsr::diff_self(x)
                  rm(list = "x"); gc(); gc();
                  mask_img(r, mask)
                },
                .progress = "text")
    mapply(function(img, fname){
      writenii(img, filename = fname)
    }, res, fnames)
    rm(list = "res"); gc();
  }

  # ret_mask_fname
  # les_out_fname
  if (verbose) {
    message("Saving List of Names")
  }
  L = list(mask = ret_mask_fname,
           all_imgs = every_fname,
           outdir = outdir)
  if (have_gold_standard) {
    L$gold_standard = les_out_fname
  }
  if (verbose) {
    message("Returning List of Names")
  }
  return(L)

}
