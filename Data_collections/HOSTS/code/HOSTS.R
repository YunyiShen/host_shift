get_HOSTS_table <- function(w){
  html_table(w[[5]])
}

library(rvest)
library(tidyverse)
base_url <- "https://www.nhm.ac.uk/our-science/data/hostplants/search/list.dsml?Genusqtype=starts+with&Familyqtype=starts+with&PFamilyqtype=starts+with&PGenusqtype=starts+with&sort=Family&beginIndex=<page>&PSpeciesqtype=starts+with&Speciesqtype=starts+with&searchPageURL=index%2edsml%3famp%3bamp%3bFamilyqtype%3dstarts%2bwith%26amp%3bamp%3bGenusqtype%3dstarts%2bwith%26Speciesqtype%3dstarts%2bwith%26amp%3bFamilyqtype%3dstarts%2bwith%26amp%3bSpeciesqtype%3dstarts%2bwith%26amp%3bamp%3bPGenusqtype%3dstarts%2bwith%26PFamilyqtype%3dstarts%2bwith%26Genusqtype%3dstarts%2bwith%26amp%3bGenusqtype%3dstarts%2bwith%26Familyqtype%3dstarts%2bwith%26PGenusqtype%3dstarts%2bwith%26amp%3bPFamilyqtype%3dstarts%2bwith%26amp%3bamp%3bsort%3dFamily%26amp%3bPGenusqtype%3dstarts%2bwith%26amp%3bamp%3bSpeciesqtype%3dstarts%2bwith%26sort%3dFamily%26beginIndex%3d60%26amp%3bbeginIndex%3d30%26amp%3bsort%3dFamily%26PSpeciesqtype%3dstarts%2bwith%26amp%3bamp%3bPSpeciesqtype%3dstarts%2bwith"
curr_url <- sub("<page>","0",base_url)

page_list <- seq(0,140460,by = 30)


tt <- read_html(curr_url) %>%
  html_node("body") %>%
  html_nodes("table") %>%
  get_HOSTS_table()
