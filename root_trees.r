library(ape)
library(optparse)

################ Setting up parser for user input######
option_list <- list(
  make_option(c("-o", "--outgroup"), type="character", default=FALSE,
              help="Give the name of the outgroup by which all trees will be rooted")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)

args <- parse_args(parser, positional_arguments = 1)
opt <- args$options
file <- args$args

########################################################

#Call the variables
path = (normalizePath(dirname(file),"/", file))
outgroup1 <- opt$outgroup

#Some screen output
sprintf('Setting workdir to %s', path)
sprintf('Using %s as outgroup', outgroup1)

setwd(path)


############## Defining function to iterate through input file ################
processFile = function(filepath) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    tree <- read.tree(text=line)
    x <- root(tree, outgroup = outgroup1, resolve.root = TRUE)
    write.tree(x, file = "1000-bootreps-rooted.tree", append = TRUE, tree.names = FALSE)
  }
  close(con)
}
################################################################################

#Applying function
print("Rooting trees ...")
processFile("1000-bootreps.tree")




#This was the previous script, which works but very slow because it reads the whole file into the memory first

#print("Reading input file ...")
#input<-read.tree("1000-bootreps.tree")
#print("Rooting trees ...")
#
#list <- list(seq(0, length(input), 1000))
#
#for(i in 1:length(input)) {
#  if (i %in% list[[1]]) {
#    sprintf('Rooting tree %s', i)
#  }
#  x <- root(input[[i]], outgroup = outgroup1, resolve.root = TRUE)
#  write.tree(x, file = "1000-bootreps-rooted.tree", append = TRUE, tree.names = FALSE)
#}

