library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
input_folder <- "../../data/extracted"
output_folder <- "../../data/processed"
file_list <- list.files(path=input_folder,pattern="csv")
F_file_list <- file_list[str_detect(file_list,"forestinvariant")]
DF_file_list <- file_list[str_detect(file_list,"deforested")]
total_data <- tibble()
for (file in F_file_list){
  a <- str_split(tools::file_path_sans_ext(file),"_")%>%unlist
  stab_period=a[3]
  stab_pol=a[4]
  tmpf <- tempfile()
  x <- readLines(paste0(input_folder,"/",file))
  y <-  substr(x,2,str_length(x)-1)  
  z <- gsub( "]", "[", y )
  cat(z, file=tmpf, sep="\n")
  #leitura
#  col_names <- c("id","LIA","coords","NDVI_after","NDVI_before","after_dt","before_dt","detection_dt",
#                 "img_name_series","orig_series","orig_f_series","stab1_series","stab2_series","stab3_series","stab4_series")
  col_names <- c("id","LIA","coords","NDVI_mean","random_dt",
                 "img_name_series","orig_series","orig_f_series","stab1_series","stab2_series","stab3_series","stab4_series")
  dado_orig <- read_csv(tmpf, quote="[",
                        col_names=col_names)
  max_img <- 300
  nvars_img=sprintf("p__%03i",seq(from=1,to=max_img))
  nvars_1=sprintf("orig__%03i",seq(from=1,to=max_img))
  nvars_1f=sprintf("origf__%03i",seq(from=1,to=max_img))
  nvars_2=sprintf("harmon__%03i",seq(from=1,to=max_img))
  nvars_3=sprintf("spatial__%03i",seq(from=1,to=max_img))
  nvars_4=sprintf("harmonf__%03i",seq(from=1,to=max_img))
  nvars_5=sprintf("spatialf__%03i",seq(from=1,to=max_img))
  nvars_dates=sprintf("d__%03i",seq(from=1,to=max_img))
  dado_format <- dado_orig%>%
    mutate(nimg=str_count(img_name_series,",")+1)
  extract_date <- function(x) substr(x,19,26)
  dado_format_img <- dado_format%>%
    separate(img_name_series,into=nvars_img,sep=", ",convert=TRUE)%>%
    separate(orig_series,into=nvars_1,sep=", ",convert=TRUE)%>%
    separate(orig_f_series,into=nvars_1f,sep=", ",convert=TRUE)%>%
    separate(stab1_series,into=nvars_2,sep=", ",convert=TRUE)%>%
    separate(stab2_series,into=nvars_3,sep=", ",convert=TRUE)%>%
    separate(stab3_series,into=nvars_4,sep=", ",convert=TRUE)%>%
    separate(stab4_series,into=nvars_5,sep=", ",convert=TRUE)%>%
    mutate_at(vars(starts_with("p__")),extract_date)%>%
    mutate_all(~str_replace(.,"None","-999"))
    #mutate_all(~na_if(.,"None"))
  
  convert_millis <- function(dmillis) as.Date.numeric(as.double(dmillis)/(24*60*60*1000),origin="1970-01-01")
  dado_format_img_tidy <- dado_format_img %>%
                 pivot_longer(contains("__"),
                 names_to = c(".value", "set"),
                 names_pattern = "(.+)__(.+)")%>%
                 filter(!is.na(p))%>%filter(orig!=-999)%>%
                 mutate(date=as.Date(p,"%Y%m%d"))%>%
                 mutate_at(vars(contains("adate")),convert_millis)%>%
                 mutate(stab_period=stab_period)%>%
                 mutate(stab_pol=stab_pol)%>%
                 select(id, stab_period,stab_pol,everything())
  print (paste(stab_pol, stab_period,nrow(dado_format_img_tidy)))
  total_data <- rbind(total_data,dado_format_img_tidy)
}
write_csv(total_data, paste(output_folder,"processed_TS_forest.csv",sep="/"))
total_data <- tibble()
for (file in DF_file_list){
  a <- str_split(tools::file_path_sans_ext(file),"_")%>%unlist
  stab_period=a[3]
  stab_pol=a[4]
  tmpf <- tempfile()
  x <- readLines(paste0(input_folder,"/",file))
  y <-  substr(x,2,str_length(x)-1)  
  z <- gsub( "]", "[", y )
  cat(z, file=tmpf, sep="\n")
  #leitura
    col_names <- c("id","LIA","coords","NDVI_after","NDVI_before","after_dt","before_dt","detection_dt",
                   "img_name_series","orig_series","orig_f_series","stab1_series","stab2_series","stab3_series","stab4_series")
  #col_names <- c("id","LIA","coords","NDVI_mean","random_dt",
  #               "img_name_series","orig_series","orig_f_series","stab1_series","stab2_series","stab3_series","stab4_series")
  dado_orig <- read_csv(tmpf, quote="[",
                        col_names=col_names)
  max_img <- 300
  nvars_img=sprintf("p__%03i",seq(from=1,to=max_img))
  nvars_1=sprintf("orig__%03i",seq(from=1,to=max_img))
  nvars_1f=sprintf("origf__%03i",seq(from=1,to=max_img))
  nvars_2=sprintf("harmon__%03i",seq(from=1,to=max_img))
  nvars_3=sprintf("spatial__%03i",seq(from=1,to=max_img))
  nvars_4=sprintf("harmonf__%03i",seq(from=1,to=max_img))
  nvars_5=sprintf("spatialf__%03i",seq(from=1,to=max_img))
  nvars_dates=sprintf("d__%03i",seq(from=1,to=max_img))
  dado_format <- dado_orig%>%
    mutate(nimg=str_count(img_name_series,",")+1)
  extract_date <- function(x) substr(x,19,26)
  dado_format_img <- dado_format%>%
    separate(img_name_series,into=nvars_img,sep=", ",convert=TRUE)%>%
    separate(orig_series,into=nvars_1,sep=", ",convert=TRUE)%>%
    separate(orig_f_series,into=nvars_1f,sep=", ",convert=TRUE)%>%
    separate(stab1_series,into=nvars_2,sep=", ",convert=TRUE)%>%
    separate(stab2_series,into=nvars_3,sep=", ",convert=TRUE)%>%
    separate(stab3_series,into=nvars_4,sep=", ",convert=TRUE)%>%
    separate(stab4_series,into=nvars_5,sep=", ",convert=TRUE)%>%
    mutate_at(vars(starts_with("p__")),extract_date)%>%
    mutate_all(~str_replace(.,"None","-999"))
  #mutate_all(~na_if(.,"None"))
  
  convert_millis <- function(dmillis) as.Date.numeric(as.double(dmillis)/(24*60*60*1000),origin="1970-01-01")
  dado_format_img_tidy <- dado_format_img %>%
    pivot_longer(contains("__"),
                 names_to = c(".value", "set"),
                 names_pattern = "(.+)__(.+)")%>%
    filter(!is.na(p))%>%filter(orig!=-999)%>%
    mutate(date=as.Date(p,"%Y%m%d"))%>%
    mutate_at(vars(contains("adate")),convert_millis)%>%
    mutate(stab_period=stab_period)%>%
    mutate(stab_pol=stab_pol)%>%
    select(id, stab_period,stab_pol,everything())
  print (paste(stab_pol, stab_period,nrow(dado_format_img_tidy)))
  total_data <- rbind(total_data,dado_format_img_tidy)
}
write_csv(total_data, paste(output_folder,"processed_TS_deforested.csv",sep="/"))
