\name{extract.lib}
\alias{extract.lib.from.zip}
\alias{extract.lib.from.directory}
\alias{extract.library.tags}
\alias{compute.unique.tags}
\alias{combine.libs}
\alias{remove.sage.artifacts}
\alias{read.phd.file}
\alias{read.seq.qual.filepair}
\alias{extract.ditags}
\alias{reestimate.lib.from.tagcounts}
\alias{SAGEartifacts}
\title{Functions for SAGE library extraction}
\description{
  Functions to extract the tags in a library from sequences or
  base-caller output.
  
}
\usage{
extract.lib.from.zip(zipfile, libname=sub(".zip","",basename(zipfile)),
                     ...) 
extract.lib.from.directory(dirname, libname=basename(dirname),
                           pattern="(\.phd\.1\$)|(\.seq\$)", ...)
extract.library.tags(filelist, base.caller.format="phd",
                     remove.duplicate.ditags=TRUE, 
                     remove.N=FALSE, remove.low.quality=10,
                     taglength=10, min.ditag.length=(2*taglength-2),
                     max.ditag.length=(2*taglength+4),
                     cut.site="catg", default.quality=NA, verbose=TRUE,
                     ...) 
reestimate.lib.from.tagcounts(tagcounts, libname, default.quality=20, ...) 
compute.unique.tags(lib)
combine.libs(..., artifacts=c("Linker", "Ribosomal", "Mitochondrial"))
remove.sage.artifacts(lib,
                      artifacts=c("Linker","Ribosomal","Mitochondrial"),
                      ...)
read.phd.file(file)
read.seq.qual.filepair(file, default.quality=NA)
extract.ditags(sequence, taglength=10, filename=NA,
               min.ditag.length=(2*taglength-2),
               max.ditag.length=(2*taglength+4), cut.site="catg")

}
\arguments{
  \item{zipfile,dirname}{Name of a ZIP file or a directory that contains
    base-caller output files}
  \item{libname}{\code{libname} a character string to be assigned as
    library name}
  \item{pattern}{Regular expression to specify pattern for the files
    that will be read}
  \item{filelist}{List of files to be read}
  \item{base.caller.format}{\code{base.caller.format} can be "phd" or
  "seq" or a character vector of the length of the filelist}
  \item{remove.duplicate.ditags}{Remove duplicate ditags. TRUE or FALSE}
  \item{remove.N}{Remove all tags that contain N. TRUE or FALSE}
  \item{remove.low.quality}{Remove all tags with an average quality
  score of less than \code{remove.low.quality}. Skipped if < 0}
  \item{taglength}{Length of tags. Usually 10 or 17}
  \item{min.ditag.length,max.ditag.length}{Minimum and maximum length
    for ditags}
  \item{cut.site}{Restriction enzyme cut site. Usually CATG}
  \item{verbose}{Display information during process}
  \item{lib}{Library object}
  \item{file,filename}{Character string indicating file name}
  \item{default.quality}{Quality value to use on sequences, if quality
    files are missing}
  \item{sequence}{Construct containing sequence and quality values
    returned by read.phd.file or read.seq.qual.filepair}
  \item{artifacts}{Types of artificially generated tags to remove.}
  \item{\dots}{Arguments passed on to extraction functions.}
  \item{tagcounts}{Tagcounts from library. Integer Vecotor with Tag sequences as names.}
}
\value{
  \code{lib} returns an SAGE library object.
}
\details{
  The functions \code{extract.lib.from.zip} or
  \code{extract.lib.from.directory} should be used to extract the SAGE
  TAGS from the sequences of a library, the sequences need to be
  provided by the output files from the base caller software either
  in a ZIP archive or in a directory. These are usually the only functions
  that should directly be called by the user. The other functions are
  called by these and should only be used directly by experienced users
  to get more direct control over the process. Most arguments are
  passed on and can be specified in the high level
  functions. Zipfilenames must be specified using relative pathnames!
  
}
\references{\url{http://tagcalling.mbgproject.org}}
\author{Tim Beissbarth}

\seealso{\code{\link[sagenhaft:sage.library]{sage.library}},
  \code{\link[sagenhaft:error.correction]{error.correction}}}

\examples{
library(sagenhaft)
file.copy(system.file("data/E15postHFI.zip",package="sagenhaft"),
          "E15postHFI.zip")
E15post<-extract.lib.from.zip("E15postHFI.zip", taglength=10,
                              min.ditag.length=20, max.ditag.length=24)
E15post
}
\keyword{manip}
