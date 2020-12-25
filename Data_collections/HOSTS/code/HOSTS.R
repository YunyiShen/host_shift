get_HOSTS_table <- function(w){
  html_table(w[[5]])
}

library(rvest)
library(tidyverse)
base_url <- "https://www.nhm.ac.uk/our-science/data/hostplants/search/list.dsml?Genusqtype=starts+with&Familyqtype=starts+with&PFamilyqtype=starts+with&PGenusqtype=starts+with&sort=Family&beginIndex=<page>&PSpeciesqtype=starts+with&Speciesqtype=starts+with&searchPageURL=index%2edsml%3famp%3bamp%3bFamilyqtype%3dstarts%2bwith%26amp%3bamp%3bGenusqtype%3dstarts%2bwith%26Speciesqtype%3dstarts%2bwith%26amp%3bFamilyqtype%3dstarts%2bwith%26amp%3bSpeciesqtype%3dstarts%2bwith%26amp%3bamp%3bPGenusqtype%3dstarts%2bwith%26PFamilyqtype%3dstarts%2bwith%26Genusqtype%3dstarts%2bwith%26amp%3bGenusqtype%3dstarts%2bwith%26Familyqtype%3dstarts%2bwith%26PGenusqtype%3dstarts%2bwith%26amp%3bPFamilyqtype%3dstarts%2bwith%26amp%3bamp%3bsort%3dFamily%26amp%3bPGenusqtype%3dstarts%2bwith%26amp%3bamp%3bSpeciesqtype%3dstarts%2bwith%26sort%3dFamily%26beginIndex%3d60%26amp%3bbeginIndex%3d30%26amp%3bsort%3dFamily%26PSpeciesqtype%3dstarts%2bwith%26amp%3bamp%3bPSpeciesqtype%3dstarts%2bwith"
temp_file_path <- "./temp/<page>.csv"
log_file_path <- "./log/log.csv"
retry <- 5
sleep <- 5

page_list <- seq(0,140460,by = 30)
#page_list <- page_list[1:5]
log_cont <- data.frame(page = 0, message = "")
log_cont <- log_cont[-1,]

for(p in page_list){
  cat("start page:",p,"\n")
  curr_url <- sub("<page>",p,base_url)
  
  for(i in 1:retry){
    temp <- tryCatch( read_html(curr_url), error = function(e){e})
    Sys.sleep(sleep)
    if(!"error" %in% class(temp)){
      cat("  read success\n")
      table_temp <- temp %>%
        html_node("body") %>%
        html_nodes("table") %>%
        get_HOSTS_table()
      
      write.csv(table_temp,sub("<page>",p,temp_file_path),row.names = F)
      break
    }
    else{
      cat(" read fail\n")
      log_cont <- rbind(log_cont,data.frame(page = p,message = temp$message))
      write.csv(log_cont,log_file_path)
    }
    cat("  retry:",i+1,"\n")
  }
}

