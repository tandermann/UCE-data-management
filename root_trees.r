library(ape)
library(optparse)

################ Setting up parser for user input######
option_list <- list(
  make_option(c("-o", "--outgroup"), type="character", default=FALSE,
              help="Give the name of the outgroup by which all trees will be rooted"),
  make_option(c("-r", "--rm_taxa"), type="character", default=FALSE,
              help="List all taxa (separated by ,) that you want to exclude before rooting the trees")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)

args <- parse_args(parser, positional_arguments = 1)
opt <- args$options
file <- args$args
########################################################

#Call the variables
path = (normalizePath(dirname(file),"/", file))
outgroup1 <- opt$outgroup
remove_taxa <- opt$rm_taxa
rm_taxa <- strsplit(remove_taxa,',')
input_file = basename(file)
output_file = "1000-bootreps-rooted.tree"

#Some screen output
sprintf('Setting workdir to %s', path)
sprintf('Using %s as outgroup', outgroup1)
sprintf('Removing taxa: %s', opt$rm_taxa)

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
    if(opt$rm_taxa != FALSE) {
      tree <- drop.tip(tree, rm_taxa[[1]], trim.internal = TRUE, subtree = FALSE)
    }
    x <- root(tree, outgroup = outgroup1, resolve.root = TRUE)
    print(x)
    write.tree(x, file = output_file, append = TRUE, tree.names = FALSE)
  }
  close(con)
}
################################################################################

#Applying function
print("Rooting trees ...")
#make sure the output file is empty before starting to write
close(file(output_file, open="w"))
#Now start the function and write the rooted trees into the output file
processFile(input_file)




#This was the previous script, which works but very slow because it reads the whole file into the memory first

#path <-"/Users/tobias/GitHub/topaza_uce/preparing_for_mpest/data/topaza_allele_bootstrap"
#outgroup1 <- "Flor_0"
#setwd(path)
#print("Reading input file ...")
#input<-read.tree("test.tree")
#print("Rooting trees ...")
#
#list <- list(seq(0, length(input), 10))
#
#remove_taxa <- "Flor_1,T_pella6_0"
#rm_taxa <- strsplit(remove_taxa,',')
#
#for(i in 1:length(input)) {
#  if (i %in% list[[1]]) {
#    sprintf('Rooting tree %s', i)
#  }
#  tree <- input[[i]]
#  tree <- drop.tip(tree, rm_taxa[[1]], trim.internal = TRUE, subtree = FALSE)
#  x <- root(tree, outgroup = outgroup1, resolve.root = TRUE)
#  write.tree(x, file = "test-rooted.tree", append = TRUE, tree.names = FALSE)
#}
#
#plot(tree)
