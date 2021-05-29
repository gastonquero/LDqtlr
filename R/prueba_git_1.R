library(devtools)
library(tidyverse)
library(fs)
library(usethis)
use_git()

dir_info(all = TRUE, regexp = "^[.]git$") %>%
  select(path, type)


load_all()
