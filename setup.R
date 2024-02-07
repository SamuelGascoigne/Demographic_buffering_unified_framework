# Set-up ========

### Environment set-up ---- 

# Clear workspace

rm(list = ls()) 



### Load packages ----

# This function loads/installs the specified package as string argument.

install_load_package <- function(x) {
  
  if (!require(x, character.only = TRUE))
    
    install.packages(x)
  
  require(x, character.only = TRUE)
  
}




# A vector of the necessary packages for this script along with their function.

package_vec <- c("ggplot2",
                 "popbio",
                 "popdemo",
                 "ggridges",
                 "viridis",
                 "reshape2",
                 "scales",
                 "demogR",
                 "tidyverse",
                 "magick")

# Use the install_load_package function to call/install all of the specified packages.

sapply(package_vec, install_load_package)
