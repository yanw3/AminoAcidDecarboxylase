
library(dplyr)
library(tableHTML)
library(stringr)
library(stringi) # function stri_replace_all_regex
library(argparse)
library(gplots) # function col2hex

# Handle command line arguments first
parser <- ArgumentParser()


parser$add_argument("--workdir", type="character", default = "" ,
    help="Set the working director")

parser$add_argument("--loadDa", type="character", default = "" ,
    help="load the data previously saved")

parser$add_argument("--contextDa", type="character", default = "" ,
    help="The path to the context dataset. The dataset should have two columns. This first it the id of the neighborhood region and the second is the context")
parser$add_argument("--word_bg", type="character", default = "" ,
    help="load the eggnod list containing eggnod ids whose background will be colored")
parser$add_argument("--word_col", type="character", default = "" ,
    help="load the eggnod list containing eggnod ids whose text will be colored")


parser$add_argument("--preColDa", type="character", default = "" ,
    help="If --preCol is set TRUE, a path to previously set color panel should be provided. The panel data should be in .Rdata format")


parser$add_argument("--eggnodID", type="character", default = "" ,
    help="A file containing a list of eggnod ID that should be included simultaneously in the context.")
parser$add_argument("--eggnodIDbool", type="character", default = "all",
    help="Select 'all' or 'any' to determin whether filter eggnod ID by 'and' or 'or'.")
parser$add_argument("--contigID", type="character", default = "" ,
    help="A file containing a list of contig ID that should be included or excluded simultaneously in the context.")
parser$add_argument("--contigIDtype", type="character", default= "exclude",
    help="Choose from 'include' or 'exclude'. Determines whether to include or exclude the contig in the contigID list")
parser$add_argument("--nrID", type="character", default = "" ,
    help="A file containing a list of neighbor region ID that should be included or excluded simultaneously in the context.")
parser$add_argument("--nrIDtype", type="character", default= "exclude",
    help="Choose from 'include' or 'exclude'. Determines whether to include or exclude the neighborhood region in the nrID list")

parser$add_argument("--title", type="character", default= "title",
    help="The title of the html output file")
parser$add_argument("--outHTML", type="character", default= "outHTML.html",
    help="output of the html file")

parser$add_argument("--SaveDa", type="character", default = "" ,
    help="A path should be gave to save the color panel. The panel data should be in .Rdata format")

args <- parser$parse_args()




setwd(args$workdir)

if(nchar(args$eggnodID)>0){
  eggnodID <- read.table(args$eggnodID,sep="\t",header=FALSE,stringsAsFactors = FALSE)
}

if(nchar(args$contigID)>0){
  contigID <- read.table(args$contigID,sep="\t",header=FALSE,stringsAsFactors = FALSE)
}

if(nchar(args$nrID)>0){
  nrIDda <- read.table(args$nrID,sep="\t",header=FALSE,stringsAsFactors = FALSE)
}
# import data #
#  ------------------------------------------------------------------------------------------------------------------------------------  #
if(nchar(args$loadDa) == 0){

# import the gene context
da <- read.table(args$contextDa,sep="\t",header=FALSE,stringsAsFactors = FALSE)
# combine id and context information
da2 <- as.data.frame(paste(da[,1],da[,2],sep=":"))
names(da2) <- "context"

# import the COG ID list
word_bg <- read.table(args$word_bg,sep="\t",header=FALSE,stringsAsFactors = FALSE)
word_col <- read.table(args$word_col,sep="\t",header=FALSE,stringsAsFactors = FALSE)


# assign color #
#  ------------------------------------------------------------------------------------------------------------------------------------  #
  color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] # maximum 433
  co <- as.data.frame(color)
  co$color <- as.character((co$color))
  co_1 <- co[!grepl("white",co$color) & !grepl("black",co$color),]
  ## change 400 to a value smaller than 430
  mycol=as.data.frame(sample(co_1, 400))  

  wordcol_bg <- as.data.frame(cbind(word_bg,mycol))
  wordcol_col <- as.data.frame(cbind(word_col,mycol))
  names(wordcol_bg) <- c("id","color")
  names(wordcol_col) <- c("id","color")

  ## change color name to hex
  wordcol_bg$colhex <- col2hex(wordcol_bg$color)
  wordcol_col$colhex <- col2hex(wordcol_col$color)

  ## change bg color for COG0076@1 as red
  if(nrow(wordcol_bg[wordcol_bg$colhex == "#FF0000",]) > 0){
    wordcol_bg[wordcol_bg$colhex == "#FF0000",]$colhex <- wordcol_bg[wordcol_bg$id == "COG0076@1",]$colhex
    wordcol_bg[wordcol_bg$id == "COG0076@1",]$colhex <- "#FF0000"
    }else{
    wordcol_bg[wordcol_bg$id == "COG0076@1",]$colhex <- "#FF0000"
  }

  # link COG ID to color to generate a list for color panel #
  word_colour_bg <- list()
  for(i in 1:nrow(wordcol_bg)){
    word_colour_bg[[i]] <- wordcol_bg[i,3]
    names(word_colour_bg)[[i]] <- wordcol_bg[i,1]
  }

  word_colour_col <- list()
  for(i in 1:nrow(wordcol_col)){
    word_colour_col[[i]] <- wordcol_col[i,3]
    names(word_colour_col)[[i]] <- wordcol_col[i,1]
  }

## change overall setting: size, family, weight
da2$context <- paste("<h1 style=\"font-size:5;font-weight:30;font-family:'Arial'\">",da2$context,"</h1>",sep="")

da2_plot <- da2

}else {
    load(file=args$loadDa)
}





# plot



# da2_plot <- as.data.frame(da2[1:100,])

names(da2_plot) <- paste0("group_",gsub(";","_",args$title))

if(nchar(args$eggnodID)>0){
  if(args$eggnodIDbool == "all"){
    da2_plot <- subset(da2_plot,apply(sapply(X = eggnodID[,1], FUN = grepl, da2_plot[,1]), MARGIN=1, FUN=all))
  }else{da2_plot <- subset(da2_plot,apply(sapply(X = eggnodID[,1], FUN = grepl, da2_plot[,1]), MARGIN=1, FUN=any))}
}

if(nchar(args$nrID)>0){
  if(args$nrIDtype == "include"){
    da2_plot <- subset(da2_plot, str_detect(da2_plot[,1],paste0(nrIDda[,1],collapse="|"),negate=FALSE))
  }else{da2_plot <- subset(da2_plot, str_detect(da2_plot[,1],paste0(nrIDda[,1],collapse="|"),negate=TRUE))}
}

if(nchar(args$contigID)>0){
  if(args$contigIDtype == "include"){
    da2_plot <- subset(da2_plot, str_detect(da2_plot[,1],paste0(contigIDda[,1],collapse="|"),negate=FALSE))
  }else{da2_plot <- subset(da2_plot, str_detect(da2_plot[,1],paste0(contigIDda[,1],collapse="|"),negate=TRUE))}
}



for(i in 1:nrow(word_bg)){
da2_plot <- data.frame(lapply(da2_plot,function(x){
    gsub(word_bg[i,1],paste("<span style=\"background-color:", word_colour_bg[[word_bg[i,1]]],";font-size:5;\">", word_bg[i,1],"</span>",sep=""),x)
})) }

for(i in 1:nrow(word_col)){
da2_plot <- data.frame(lapply(da2_plot,function(x){
    gsub(word_col[i,1],paste("<span style=\"color:", word_colour_col[[word_col[i,1]]],";font-size:5;font-weight:900\">", word_col[i,1],"</span>",sep=""),x)
})) }

da2_plot %>% tableHTML(rownames = FALSE,escape = FALSE,border=0) %>% write_tableHTML(., file = args$outHTML)

# save the important data in the session to be used in the future analysis
if(nchar(args$SaveDa) > 0){
      save(da, da2, da2_plot, word_bg, word_col, word_colour_bg, word_colour_col, file=args$SaveDa)
    }
