
#' @importFrom dynwrap test_docker_installation
#' @importFrom babelwhale list_docker_images pull_container
#' @importFrom tidyr unite
#' @importFrom dplyr pull

execute_container <- function(ref_data, verbose) {

  # Process file path-----------------------------------------------------------
  ## 1. Create a temp directory to store .rds file or .h5 file for Python
  temp_dir <- tempdir()
  ## 2. Set a path for the input data
  temp_input_path <- file.path(temp_dir, "input.rds")

  # Process datasets------------------------------------------------------------
  ## 1. Convert ref_data into .rds or .h5 file and save to input path
  # simutils::write_h5files(data = ref_data, file_path = temp_input_path)
  saveRDS(ref_data, file = temp_input_path)

  # Prepare the input parameters-----------------------------------------------
  ## 1. docker image working directory
  wd <- "/home/admin/"
  ## 2. local directory of the mount point
  local_path <- temp_dir
  ## 3. docker image directory of the mount point
  docker_path <- "/home/rstudio"
  ## 4. verbose
  verbose <- verbose
  ## 5. args
  args <- NULL
  ## 6. command
  command <- NULL
  ## 7. container id
  ### (1. Check docker installation
  docker_install <- dynwrap::test_docker_installation()
  if(!docker_install){
    stop("Docker has not been installed or started! Please check it!")
  }
  ### (2. Check the installation of docker image
  images <- babelwhale::list_docker_images() %>%
    tidyr::unite("Repository", "Tag", sep = ":", col = "Image") %>%
    dplyr::pull("Image")

  if(!"xinern/sctreemerge:v1.0.0" %in% images){
    # If not, pull duohongrui/simpipe:latest
    babelwhale::pull_container(container_id = "xinern/sctreemerge:v1.0.0")
  }
  ### (3. docker container id
  container_id <- "xinern/sctreemerge:v1.0.0"




  # Run container---------------------------------------------------------------
  output <- babelwhale::run(container_id = container_id,
                            command = command,
                            args = args,
                            volumes = paste0(local_path, ":", docker_path),
                            workspace = wd,
                            verbose = verbose,
                            debug = FALSE)

  # Get result------------------------------------------------------------------
  if(verbose){
    cat("Output is saved to ", local_path, "\n")
    cat("Attempting to read output into R\n")
  }

  ## return output
  readRDS(file = file.path(local_path, "output.rds"))

}
