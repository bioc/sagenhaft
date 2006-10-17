
## Extraction of tags of a SAGE library

extract.lib.from.zip <- function(zipfile, libname=sub(".zip","",basename(zipfile)), ...) {
  if(.Platform$OS.type == "windows") {
    tmpdir <- tempfile()
    dir.create(tmpdir)
    zip.unpack(zipfile, tmpdir)
    lib <- extract.lib.from.directory(tmpdir, libname, ...)
    unlink(tmpdir, recursive=TRUE)
  } else {
    unlink(tempdir(), recursive=TRUE)
    dir.create(tempdir())
    tmpdir <- zip.file.extract("", zipfile)
    lib <- extract.lib.from.directory(tmpdir, libname, ...)
    unlink(tmpdir, recursive=TRUE)
    dir.create(tempdir())
  }
  return(lib)
}

extract.lib.from.directory <- function(dirname, libname=basename(dirname), pattern="(\\.phd\\.1$)|(\\.seq$)", ...) {
  filelist <- dir(dirname, pattern, full.names=TRUE)  
  if(length(filelist)==0) { stop("No files!") }
  base.caller.format <- rep("phd", length(filelist))
  base.caller.format[grep("\\.seq$", filelist)] <- "seq"
  lib <- extract.library.tags(filelist, base.caller.format=base.caller.format, ...)
  if(!nrow(lib$seqs)>0) stop("No tags found in sequence files!")
  lib <- compute.unique.tags(lib)
  if(!nrow(lib$tags)>0) stop("Error in finding unique tags!")
  lib$libname <- libname
  lib <- estimate.errors.mean(lib)
  lib <- em.estimate.error.given(lib, ...)
  lib <- remove.sage.artifacts(lib, ...)
  return(lib)
}

reestimate.lib.from.tagcounts <- function(tagcounts, libname, default.quality=20, ...) {
  if(is.character(tagcounts)) { tmp <- tagcounts; tagcounts <- integer(length(tagcounts))+1; names(tagcounts)<-tmp }
  if(!is.numeric(tagcounts)) stop("Tagcounts needs to be a numeric vector with tagsequences as names or character vector of tagsequences in library!")
  if(!is.character(names(tagcounts))) stop("Tagcounts needs to have names giving the tag sequences!")
  taglength = unique(nchar(names(tagcounts)))
  if(length(taglength)!=1) stop("All tag sequences need to be of same lengths!")
  tagnames <- names(tagcounts)
  tagcounts <- as.integer(tagcounts)
  tags <- rep.int(tagnames, tagcounts)
  if(!length(tags)>0) stop("No tags found!")
  error.scores <- matrix(default.quality, nrow=length(tags), ncol=taglength)
  colnames(error.scores) <- paste("q", 1:taglength, sep="")
  seqs <- data.frame(seq=I(tags), seqextra=I(sample(c("a","c","g","t"),length(tags),TRUE)), error.scores,
                     ditaglength=rep(NA,length(tags)), file=I(rep("",length(tags))))
  comment <- c(paste("# date:", date(), sep=" "),
               "# base.calling.method: Unknown",
               paste("# default.quality:", default.quality, sep=" "),
               )
  lib <- list(libname, taglength=taglength, seqs=seqs, comment=comment)
  class(lib) <- "sage.library"
  lib <- compute.unique.tags(lib)
  if(!nrow(lib$tags)>0) stop("Error in finding unique tags!")
  lib <- estimate.errors.mean(lib)
  lib <- em.estimate.error.given(lib, ...)
  lib <- remove.sage.artifacts(lib, ...)
  return(lib)
}

combine.libs<-function(..., artifacts=c("Linker", "Ribosomal", "Mitochondrial")) {
  arglist <- list(...)
  if(length(arglist) < 2) stop("Less than two libraries to combine!")
  l <- arglist[[1]]
  if(class(l) != "sage.library") stop("Arguments must be of class sage.library!")
  ln <- list(libname=l$libname, taglength=l$taglength, seqs=l$seqs,
             comment=c(paste("# Combined Library 1: ", l$libname, sep=""), l$comment))
  class(ln) <- "sage.library"
  for(i in 2:length(arglist)) {
    l<-arglist[[i]]
    if(class(l) != "sage.library") next
    if(l$taglength != ln$taglength) stop("Libraries of different tag length cannot be combined!")
    ln$libname <- paste(ln$libname, l$libname, sep=".")
    ln$comment <- c(ln$comment, paste("# Combined Library ", i, ": ", l$libname, sep=""), l$comment)
    ln$seqs <- rbind(ln$seqs, l$seqs)
  }
  ln <- compute.unique.tags(ln)
  ln <- estimate.errors.mean(ln)
  ln <- em.estimate.error.given(ln, maxstep=50)
  ln <- remove.sage.artifacts(ln, artifacts=artifacts)
  return(ln)
}

compute.unique.tags <- function(lib) {
  tagnums <- tagsequence2tagnum(lib$seqs[,"seq"], lib$taglength)
  counts <- table(tagnums)
  tags  <- tagnum2tagsequence(as.numeric(names(counts)), lib$taglength)
  counts <- as.integer(counts)
  extra <- as.matrix(round(100*table(tagnums, as.character(lib$seqs[,"seqextra"]))/counts))
  extra <- extra[, c("a","c","g","t")]
  rownames(extra) <- NULL
  class(extra) <- NULL
  ditaglength <- as.integer(round(tapply(as.numeric(lib$seqs[,"ditaglength"]), tagnums, mean)))
  error <- as.integer(round(tapply(as.numeric(rowMeans(lib$seqs[,paste("q", 1:lib$taglength, sep="")])), tagnums, mean)))
  lib$tags <- data.frame(tag=I(tags), count.raw=counts, extra, avg.ditaglength=ditaglength, avg.error.score=error)
  lib$ntag <- length(counts)
  lib$nseq <- sum(counts)
  return(lib)
}

remove.sage.artifacts <- function(lib, artifacts=c("Linker", "Ribosomal", "Mitochondrial"), ...) {
  data(SAGEartifacts, envir=environment())
  if(lib$taglength == 10) { taglength <- "Short" } else { taglength <- "Long" }
  for(a in artifacts) {
    artifact.tags <- tolower(as.character(SAGEartifacts[as.character(SAGEartifacts[,2])==paste(taglength, a, sep=""), 1]))
    artifact.tags <- match(artifact.tags, lib$tags[,"tag"])
    artifact.tags <- artifact.tags[!is.na(artifact.tags)]
    if(length(artifact.tags)<1) next
    lib$comment <- c(lib$comment, paste("# Removed ", taglength, a, ": ", sum(lib$tags[artifact.tags,"count.raw"]), " ",
                                        sum(lib$tags[artifact.tags,"count.adjusted"]), sep=""))
    lib$tags <- lib$tags[-artifact.tags,]
  }
  return(lib)
}

extract.library.tags <- function(filelist, base.caller.format="phd", remove.duplicate.ditags=TRUE, remove.N=FALSE, remove.low.quality=10,
                                 taglength=10, min.ditag.length=(2*taglength-2), max.ditag.length=(2*taglength+4), cut.site="catg",
                                 default.quality=NA, verbose=TRUE, ...) {
  # remove.low.quality if not -1 remove tags where the average of the quality scores is less than remove.low.quality
  base.caller.format <- rep(base.caller.format, length.out=length(filelist))
  libditags <- NULL
  for(i in 1:length(filelist)) {
    file <- filelist[i]
    if(verbose) print(paste(i, file, sep=":"))
    if(base.caller.format[i] == "phd") {
      seq <- read.phd.file(file)
    } else if(base.caller.format[i] == "seq") {
      seq <- read.seq.qual.filepair(file, default.quality)
    } else {
      stop("file format not supported")
    }
    ditags <- extract.ditags(seq, taglength, basename(file), min.ditag.length, max.ditag.length, cut.site)
    if(!is.null(ditags) && nrow(ditags) > 0) {
      if(!is.null(libditags)) {
        libditags <- rbind(libditags, ditags)
      } else {
        libditags <- ditags
      }
    }
  }
  if(is.null(libditags)) { stop("no tags found") }
  taglength <- nchar(libditags[1,3])
  nduplicate.ditags <- nrow(libditags) - length(unique(libditags[,1]))
  if(remove.duplicate.ditags) {
    o<-order(rowSums(libditags[,c(5:(taglength+4),(taglength+7):(2*taglength+6))]), decreasing=TRUE)
    libditags <- libditags[o,]
    libditags <- libditags[match(unique(libditags[,1]), libditags[,1]),]
  }
  tags <- libditags[,c(3:(taglength+4),2, ncol(libditags))]
  tags[(nrow(tags)+1):(2*nrow(tags)),] <- libditags[,c((taglength+5):(2*taglength+6),2,ncol(libditags))]
  names(tags)[1:2] <- c("seq", "seqextra")
  names(tags)[3:(taglength+2)] <- paste("q", 1:taglength, sep="")
  if(remove.N) {
    tags <- tags[grep("^[acgt]*$", tags[,1]),]
  }
  if(remove.low.quality>0) {
    tags <- tags[rowMeans(tags[,paste("q", 1:taglength, sep="")]) >= remove.low.quality,]
  }
  rownames(tags) <- 1:nrow(tags)
  comment <- c(paste("# nfiles:", length(filelist), sep=" "),
               paste("# date:", date(), sep=" "),
               paste(c("# base.calling.method:", sort(unique(base.caller.format))), collapse=" "),
               paste("# remove.duplicate.ditags:", remove.duplicate.ditags, sep=" "),
               paste("# nduplicate.ditags:", nduplicate.ditags, sep=" "),
               paste("# remove.N:", remove.N, sep=" "),
               paste("# remove.low.quality:", remove.low.quality, sep=" "),
               paste("# min.ditag.length:", min.ditag.length, sep=" "),
               paste("# max.ditag.length:", max.ditag.length, sep=" "),
               paste("# cut.site:", cut.site, sep=" "),
               )
  lib <- list(taglength=taglength, seqs=tags, comment=comment)
  class(lib) <- "sage.library"
  return(lib)
}

read.phd.file <- function(file) {
  lines <- readLines(file)
  dna<-lines[(match("BEGIN_DNA", lines)+1):(match("END_DNA", lines)-1)]
  dna<-t(matrix(unlist(strsplit(dna, " ", FALSE)), nrow=3))
  return(list(dna=tolower(dna[,1]), error.scores=as.integer(dna[,2])))
}

read.seq.qual.filepair <- function(file, default.quality=NA) { 
  lines <- readLines(sub("\\.[^\.]*$", ".seq", file))
  if(length(lines) < 1) { stop(paste("File missing or wrong format:", sub("\\.[^\.]*$", ".seq", file))); }
  dna<-unlist(strsplit(gsub("[^A-Za-z]","", paste(lines[-1], collapse="")), ""))
  qfile <- sub("\\.[^\.]*$", ".qual", file)
  if(file.access(qfile, 4)==-1) {
    if(is.na(default.quality)) {
      stop(paste("File unreadable or missing:", qfile));
    } else {
      qual<-as.integer(rep(default.quality, length(dna)))
    }
  } else {
    lines <- readLines(qfile)
    if(length(lines) < 1) { stop(paste("File has wrong format:", qfile)); }
    qual<-as.integer(unlist(strsplit(gsub("[^0-9][^0-9]*"," ", paste(lines[-1], collapse=" ")), " ")))
  }
  if(length(dna)!=length(qual)) { stop(paste("Files have wrong format:", file)); }
  return(list(dna=tolower(dna), error.scores=qual))
}

extract.ditags <-function(sequence, taglength=10, filename=NA, min.ditag.length=(2*taglength-2), max.ditag.length=(2*taglength+4),
                          cut.site="catg") {
  s<-unlist(strsplit(gsub("catg", "CATG", paste(sequence$dna, collapse="")),""))
  n <- cumsum(s=="G")
  n[s=="C" | s=="A" | s=="T" | s=="G"] <- 0
  ditags<-split(sequence$dna, n)
  ditags <- ditags[c(-1,-length(ditags))]
  ditags <- sapply(ditags, paste, collapse="")
  error.scores <-split(sequence$error.scores,n)
  error.scores <- error.scores[c(-1,-length(error.scores))]
  ditaglength <- nchar(ditags)
  ditags <- ditags[ditaglength>=min.ditag.length & ditaglength<=max.ditag.length]
  error.scores <- error.scores[ditaglength>=min.ditag.length & ditaglength<=max.ditag.length]
  ditaglength <- ditaglength[ditaglength>=min.ditag.length & ditaglength<=max.ditag.length]
  if(length(ditaglength)==0) { return(NULL) }
  tag1 <- substr(ditags, 1, taglength)
  tag1extra <- substr(ditags, taglength+1, taglength+1)
  tag2 <- revcomp(substr(ditags, ditaglength-taglength+1, ditaglength))
  tag2extra <- revcomp(substr(ditags, ditaglength-taglength, ditaglength-taglength))
  error.scores1 <- t(sapply(error.scores, .subset, 1:taglength))
  error.scores2 <- lapply(error.scores, rev)
  error.scores2 <- t(sapply(error.scores2, .subset, 1:taglength))
  # PLAN rename score columns
  # PLAN if 2*taglength > ditaglength set quality of last base to low
  return(data.frame(ditags=I(ditags), ditaglength=ditaglength, tag1=I(tag1), tag1extra=I(tag1extra), error.scores1,
              tag2=I(tag2), tag2extra=I(tag2extra), error.scores2, check.names=FALSE, file=I(rep(filename, length(ditags)))))
}

# sequence error correction

estimate.errors.mean <- function(lib) {
  taglength <- lib$taglength
  tagnames <- tagsequence2tagnum(lib$tags[,"tag"], taglength)
  # compute 3*taglength+1 neighbors (tags with one base change)
  rtags <- compute.sequence.neighbors(lib$seqs[,"seq"], taglength, lib$seqs[,paste("q", 1:taglength, sep="")], "numeric")
  rtags.a <- rtags$quality
  rtags <- rtags$tags
  # convert tag numbers to indexes in array of tagnames                           
  # omit all tags that do no appear in the library
  unique.tag.count <- nrow(rtags)
  rtags <- match(as.vector(rtags), tagnames, NA)
  # correct errors with proportions
  rtags <- matrix(rtags, nrow=unique.tag.count)
  proportions <- lib$tags[,"count.raw"] / lib$nseq
  k<-1
  for(j in 1:taglength) {
    proportions.sum <- rowSums(cbind(0, proportions[rtags[,k+1]], proportions[rtags[,k+2]], proportions[rtags[,k+3]]), na.rm=TRUE)
    rtags.a[,k+1] <- pmax(0, rtags.a[,k+1] * proportions[rtags[,k+1]] / proportions.sum, na.rm=TRUE)
    rtags.a[,k+2] <- pmax(0, rtags.a[,k+2] * proportions[rtags[,k+2]] / proportions.sum, na.rm=TRUE)
    rtags.a[,k+3] <- pmax(0, rtags.a[,k+3] * proportions[rtags[,k+3]] / proportions.sum, na.rm=TRUE)
    k <- k + 3
  }
  rtags.a[,1] <- pmax(0.01, pmin(1, 1 - rowSums(rtags.a[,-1], na.rm=TRUE), na.rm=TRUE))
  rtags.a <- rtags.a/rowSums(rtags.a)     # normalize row sums to 1 ???
  
  # use sparse array
  unique.all.tag.count <- length(tagnames)
  lib$alpha <- create.matrix.csr(as.numeric(rtags.a),
                                 rep(rtags[1:unique.tag.count], 1+taglength*3), rtags,
                                 dim=c(unique.all.tag.count, unique.all.tag.count))
  return(lib)
}

compute.sequence.neighbors <- function(tags, taglength=10, quality.scores=NULL, output="character") {
  if(is.character(tags)) {
    tags <- tagsequence2tagnum(tags, taglength)
  } else {
    tags <- as.numeric(tags)
  }
  names(tags) <- NULL
  unique.tags <- sort(as.numeric(unique(tags)))
  unique.tag.count <- length(unique.tags)  
  quality.scores <- 10^(-as.matrix(quality.scores)/10)
  rownames(quality.scores) <- NULL
  colnames(quality.scores) <- NULL
  
  # convert numbers array representing tags
  tags.matrix <- tagnum2tagmatrix(unique.tags, taglength)
  rtags <- matrix(0, nrow=unique.tag.count, ncol=1+taglength*3)
  rtags.quality <- rtags
  k<-1
  prob <- rep(NA, unique.tag.count)

  # compute 3*taglength+1 neighbors (tags with one base change)
  for(j in 1:taglength) {
    # take means for probabilities
    if(!is.null(quality.scores)) prob <- tapply(quality.scores[,j], tags, mean)
    prob <- prob/3
    for(b in 1:3) {
      k<-k+1
      # do one base exchange in position j, bases are 0:3 (acgt),
      # here if b is the original base it is replaced with 0, otherwise with b
      tags.matrix.tmp <- tags.matrix
      temp <- b==tags.matrix[,j]
      tags.matrix.tmp[temp,j] <- 0
      tags.matrix.tmp[!temp,j] <- b
      # convert tags to numbers
      rtags[,k]<-tagmatrix2tagnum(tags.matrix.tmp, taglength)
      # assign probabiliies
      rtags.quality[,k] <- prob
    }
  }
  rtags[,1] <- unique.tags
  rtags.quality[,1] <- pmax(0.01, 1 - rowSums(rtags.quality[,-1]))
  rtags.quality <- rtags.quality/rowSums(rtags.quality)     # normalize row sums to 1 ???
  if(output=="character") rtags <- matrix(tagnum2tagsequence(rtags, taglength), nrow=unique.tag.count)
  return(list(tags=rtags, quality=rtags.quality))
}

em.estimate.error.given<-function(lib, maxstep=50, ...) {

  library(SparseM)
  # initial complete data
  lambda <- sum(lib$nseq)
  m <- lib$tags[,"count.raw"]
  n <- m
  l <- numeric(maxstep)
  v <- numeric(maxstep)
  at <- t(lib$alpha)
  
  for(i in 1:maxstep) {
    pm <- m

    # estimate new parameters
    p <- m/lambda

    # compute E values
    rsums <- as.vector(n / (lib$alpha %*% p))
    rsums[is.na(rsums) | is.infinite(rsums)] <- 0
    m <- p * as.vector(at %*% rsums)
    
    # compute likelihood
    l[i] <- -lambda + sum(m*log(p*lambda), na.rm=TRUE)

    # compute difference to real
    if("count.true" %in% names(lib$tags)) {
      v[i] <- sum(abs(m-lib$tags[,"count.true"]))
    } else {
      v[i] <- sum(abs(m-n))
    }
    
    if(all(pm==m)) {
      break
    }
  }
  if(i>10) {
    lib$comment <- c(lib$comment,
                     paste("# EM steps:", i, sep=" "),
                     paste(c("# likelihood (every 10 steps):", round(l[c(1,seq(10,i-1,10),i)],1)), collapse=" "),
                     paste(c("# var (every 10 steps):", round(v[c(1,seq(10,i-1,10),i)],1)), collapse=" "))
  }
  lib$likelihood <- l[1:i]
  lib$var <- v[1:i]
  names(m) <- NULL
  names(p) <- NULL
  lib$tags$count.adjusted <- round(m,2)
  lib$tags$prop.estimate <- p

  return(lib); 
}

# simulation of SAGE libraries

sagelibrary.simulate<-function(taglength=4, lambda=1000, mean.error=0.01, error.sd=1, withintagerror.sd=0.2,
                      ngenes=min(4^taglength, 100000), base.lib=NULL, libseed=-1, ...) {
  #maybe better mean.real=0.0002559957 median.real=0.0001584893 sd.accross=2.288200 sd.within=0.9313946

  ntags <- 4^taglength
  if(is.null(base.lib)) {
    if(libseed==-1) libseed<-10*floor(runif(1, 1, 10000))
    libname <- sprintf("SIM%5.5dL%2.2dE%2.2dC%d", as.integer(libseed), as.integer(taglength),
                       as.integer(round(-10*log10(mean.error))), as.integer(lambda))
    set.seed(libseed)
    if(ngenes < ntags) {
      genenames <- unique(ceiling(runif(ngenes, 0, ntags)))
    } else {
      genenames <- 1:ntags      
    }
    p <-  rlnorm(ngenes, log(0.1), 2)
    # maybe better real.mean=8.8e-06 real.median=5.3e-06 real.sd=0.83
    p <- p / sum(p)
    diff.true <- NULL
  } else {
    # construct library based on other simulated library with some known differences
    if(libseed==-1) libseed<-as.integer(sub("^SIM([0-9]*).*$", "\\1", base.lib$libname))+1
    libname <- sprintf("SIM%5.5dL%2.2dE%2.2dC%d", as.integer(libseed), as.integer(taglength),
                       as.integer(round(-10*log10(mean.error))), as.integer(lambda))
    set.seed(libseed)
    genenames <- base.lib$all.genes
    p <- base.lib$prop.true
    diff.true <- rep(c(1,1,1,1,1,2,-1,1,1,1,0.5,-1,1,1,1,3,-2,-3,1,1,1,0.33333,-1,1,1,1,1), length.out=length(p))    
    o<- order(p)
    pn <- p[o]*diff.true
    pn[diff.true==-1] <- (p[o])[diff.true==-1] - pn[which(diff.true==-1)-1] + (p[o])[which(diff.true==-1)-1]
    pn[diff.true==-2] <- (p[o])[diff.true==-2] - (pn[which(diff.true==-2)-1] - (p[o])[which(diff.true==-2)-1])/2
    pn[diff.true==-3] <- (p[o])[diff.true==-3] - (pn[which(diff.true==-3)-2] - (p[o])[which(diff.true==-3)-2])/2
    pn <- pn[order(o)]/sum(pn)    
    diff.true <- pn / p
    p <- pn    
  }
  ngenes<-length(genenames)
  m <- as.integer(rpois(ngenes, p*lambda))
  tagnames<-which(m>0)
  po<- p[tagnames]
  m<- m[tagnames]
  tagnames <- genenames[tagnames]
  tagcount <- sum(m)

  # create matrix with random errors
  # tag sequence from number
  a.real <- rep(tagnames, m)
  a.real.matrix <- tagnum2tagmatrix(a.real, taglength)
  # generate errors
  a <- matrix(rlnorm(taglength*tagcount, log(rlnorm(tagcount, log(mean.error), error.sd)),withintagerror.sd), ncol=taglength)
  a[a>1]<-1
  rtags <- matrix(0, nrow=tagcount, ncol=taglength)
  for(j in 1:taglength) {
    prob<-cbind(1-a[,j], a[,j]/3, a[,j]/3, a[,j]/3)
    # sample bases
    rtags[,j] <- apply(prob, 1, sample, x=0:3, size=1, replace=TRUE)
  }
  swap <- a.real.matrix == rtags
  rtags[rtags==0] <- a.real.matrix[rtags==0]
  rtags[swap] <- 0
  # tag sequences to numbers
  a.names <- tagmatrix2tagnum(rtags, taglength)
  n <- table(a.names)
  tagnames <- unique(c(tagnames, as.numeric(names(n))))
  if(length(tagnames) > length(m)) {
    m[(length(m)+1):length(tagnames)] <- 0
    po <- p[match(tagnames, genenames, NA)]
  }
  n <-as.integer(n[match(tagnames, as.numeric(names(n)), NA)])
  n[is.na(n)] <- 0

  # Use Sparse arrays
  unique.tag.count <- length(tagnames)
  net <- table.sparse(a.real, a.names)
  ne <- create.matrix.csr(net$count, match(net$a.names, tagnames, NA), match(net$a.real, tagnames, NA),
                          dim=c(unique.tag.count,unique.tag.count))
  alpha<-ne
  alpha@ra <- alpha@ra / m[alpha@ja]
  # format in sage.library format
  a <- round(-10 * log10(a))
  error <- round(tapply(rowMeans(a), a.names, mean))
  error <- as.integer(error[match(as.numeric(tagnames),as.numeric(names(error)),NA)])
  colnames(a) <- paste("q", 1:ncol(a), sep="")
  rownames(a) <- NULL  
  tags <- data.frame(tag=I(tagnum2tagsequence(tagnames, taglength)), count.raw=n, count.true=m, prop.true=po, avg.error.score=error)
  seqs <- data.frame(seq=I(tagnum2tagsequence(a.names, taglength)), seq.true=I(tagnum2tagsequence(a.real, taglength)))
  seqs <- cbind(seqs, a)
  if(!is.null(diff.true)) {
     diff.true <- as.numeric(diff.true[match(tagnames, genenames, NA)])
     tags <- cbind(tags, diff.true)
  }
  comment <- c(paste("# lambda:", lambda, sep=" "),
               paste("# date:", date(), sep=" "),
               paste("# ngenes:", ngenes, sep=" "),
               paste("# seed:", libseed, sep=" "),
               paste("# mean.error:", mean.error, sep=" "),
               paste("# error.sd:", error.sd, sep=" "),
               paste("# withintagerror.sd:", withintagerror.sd, sep=" "),
               )
  lib <- list(libname=libname, taglength=taglength, ntags=unique.tag.count, nseq=tagcount, tags=tags, seqs=seqs,
              count.table=ne, alpha.true=alpha, comment=comment, all.genes=genenames, prop.true=p, seed=libseed)
  class(lib) <- "sage.library"
  lib <- estimate.errors.mean(lib)
  lib <- em.estimate.error.given(lib, ...)
  return(lib)
}


# READING and WRITING SAGE libraries

write.sage.library <- function(x, file=paste(x$libname, "sage", sep="."), what="complete") {
  write(c(paste("# libname:", x$libname, sep=" "),
          paste("# nseq:", x$nseq, sep=" "),
          paste("# ntag:", x$ntag, sep=" "),
          paste("# taglength:", x$taglength, sep=" "), 
          x$comment, "# library counts",
          paste(names(x$tags), collapse="\t")), file)
  write.table(x$tags, file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="", append=TRUE)
  if(what=="complete") {
    write(c("# library sequences",
            paste(names(x$seqs), collapse="\t")), file, append=TRUE)
    write.table(x$seqs, file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="", append=TRUE)
  }
}

read.sage.library <- function(file) {
  lines <- readLines(file)
  s1 <- match("# library counts", lines)
  s2 <- match("# library sequences", lines)
  s3 <- s2
  if(is.na(s1)) stop("Not a valid library file, tag counts missing!")
  if(is.na(s2)) s3 <- length(lines+1)
  comment <- lines[1:(s1-1)]
  temp <- charmatch("# libname:", comment)
  if(is.na(temp)) stop("Not a valid library file, libname missing!")
  libname <- sub("^# libname: *", "", comment[temp])
  comment <- comment[-temp]
  temp <- charmatch("# ntag:", comment)
  if(is.na(temp)) stop("Not a valid library file, ntag missing!")
  ntag <- as.integer(sub("^# ntag: *", "", comment[temp]))
  comment <- comment[-temp]
  temp <- charmatch("# nseq:", comment)
  if(is.na(temp)) stop("Not a valid library file, nseq missing!")
  nseq <- as.integer(sub("^# nseq: *", "", comment[temp]))
  comment <- comment[-temp]
  temp <- charmatch("# taglength:", comment)
  if(is.na(temp)) stop("Not a valid library file, taglength missing!")
  taglength <- as.integer(sub("^# taglength: *", "", comment[temp]))
  comment <- comment[-temp]
  tags<-read.table(file, skip=s1, nrows=s2-s1-2, header=TRUE, sep="\t", as.is=TRUE, comment.char="")
  if(is.na(s2)) {
    seqs<-NULL
  } else {
    seqs<-read.table(file, skip=s2, header=TRUE, sep="\t", as.is=TRUE, comment.char="")
  }
  lib <- list(libname=libname, nseq=nseq, ntag=ntag, taglength=taglength, tags=tags, seqs=seqs, comment=comment)
  class(lib) <- "sage.library"
  return(lib)
}

write.sage.library.comparison <- function(x, file=paste(x$name, "sagecomp", sep=".")) {
  write(c(paste("# name:", x$name, sep=" "),
          paste("# ntag:", x$ntag, sep=" "),
          paste("# taglength:", x$taglength, sep=" "), 
          x$comment, "# comparison results",
          paste(names(x$data), collapse="\t")), file)
  write.table(x$data, file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", na="", append=TRUE)
}

read.sage.library.comparison <- function(file) {
  lines <- readLines(file)
  s1 <- match("# comparison results", lines)
  if(is.na(s1)) stop("Not a valid library file, header line missing!")
  comment <- lines[1:(s1-1)]
  temp <- charmatch("# name:", comment)
  if(is.na(temp)) stop("Not a valid library file, name missing!")
  libname <- sub("^# name: *", "", comment[temp])
  comment <- comment[-temp]
  temp <- charmatch("# ntag:", comment)
  if(is.na(temp)) stop("Not a valid library file, ntag missing!")
  ntag <- as.integer(sub("^# ntag: *", "", comment[temp]))
  comment <- comment[-temp]
  temp <- charmatch("# taglength:", comment)
  if(is.na(temp)) stop("Not a valid library file, taglength missing!")
  taglength <- as.integer(sub("^# taglength: *", "", comment[temp]))
  comment <- comment[-temp]
  tags<-read.table(file, skip=s1, header=TRUE, sep="\t", as.is=TRUE, comment.char="")
  lib <- list(name=libname, ntag=ntag, taglength=taglength, data=tags, comment=comment)
  class(lib) <- "sage.library.comparison"
  return(lib)
}

# Library comparison

compare.lib.pair <- function(lib1, lib2) {
  if(lib1$taglength != lib2$taglength) stop("Libraries of different tag length not comparable!");
  taglength <- lib1$taglength
  tagnames <- sort(unique(c(lib1$tags[,"tag"], lib2$tags[,"tag"])))
  libname <- paste(lib1$libname, lib2$libname, sep=":")
  
  # match tags of library 1
  o<-match(tagnames, lib1$tags[,"tag"], NA)
  tagcount1 <- lib1$nseq
  em1 <- lib1$tags[o,"count.adjusted"]
  em1[is.na(em1)] <-0
  n1 <- lib1$tags[o,"count.raw"]
  n1[is.na(n1)] <-0
  if("count.true" %in% names(lib1$tags)) {
    m1 <- lib1$tags[o,"count.true"]
    m1[is.na(m1)] <-0
  } else {
    m1 <- NULL
  }

  # match tags of library 2
  o<-match(tagnames, lib2$tags[,"tag"], NA)
  tagcount2 <- lib2$nseq
  em2 <- lib2$tags[o,"count.adjusted"]
  em2[is.na(em2)] <-0
  n2 <- lib2$tags[o,"count.raw"]
  n2[is.na(n2)] <-0
  if("count.true" %in% names(lib2$tags)) {
    m2 <- lib2$tags[o,"count.true"]
    m2[is.na(m2)] <-0
  } else {
    m2 <- NULL
  }

  M.adjusted <- log2(((em1+0.5) * (tagcount2-em2+0.5)) / ((em2+0.5) * (tagcount1-em1+0.5)))
  tests.adjusted <- sage.test(round(em1),round(em2),n1=tagcount1, n2=tagcount2)
  M.raw <- log2(((n1+0.5) * (tagcount2-n2+0.5)) / ((n2+0.5) * (tagcount1-n1+0.5)))
  tests.raw <- sage.test(n1,n2,n1=tagcount1, n2=tagcount2)
  
  if(!is.null(m1) && !is.null(m2)) {
    A <- (log2(0.5+round(m1*(tagcount1+tagcount2)/(2*tagcount1)))+log2(0.5+round(m2*(tagcount1+tagcount2)/(2*tagcount2))))/2
    M.true <- log2(((m1+0.5) * (tagcount2-m2+0.5)) / ((m2+0.5) * (tagcount1-m1+0.5)))
    tests.true <- sage.test(m1,m2,n1=tagcount1, n2=tagcount2)
    tags <- data.frame(tag=I(tagnames), A.true=A, count1.raw=n1, count2.raw=n2, M.raw=M.raw, tests.raw=tests.raw,
                       count1.adjusted=em1, count2.adjusted=em2, M.adjusted=M.adjusted, tests.adjusted=tests.adjusted,
                       count1.true=m1, count2.true=m2, M.true=M.true, tests.true=tests.true)
  } else {
    A <- (log2(0.5+round(em1*(tagcount1+tagcount2)/(2*tagcount1)))+log2(0.5+round(em2*(tagcount1+tagcount2)/(2*tagcount2))))/2
    tags <- data.frame(tag=I(tagnames), A.adjusted=A, count1.raw=n1, count2.raw=n2, M.raw=M.raw, tests.raw=tests.raw,
                       count1.adjusted=em1, count2.adjusted=em2, M.adjusted=M.adjusted, tests.adjusted=tests.adjusted)
  }

  if("diff.true" %in% names(lib2$tags)) {
    diff.true<-lib2$tags[o,"diff.true"]
    diff.col <- numeric(length(diff.true)) + 2
    diff.col[abs(diff.true-1)<0.1] <- 1
    diff.col[abs(diff.true-2)<0.1|abs(diff.true-0.5)<0.01] <- 3
    diff.col[abs(diff.true-3)<0.1|abs(diff.true-0.33333)<0.01] <- 4
    tags <- cbind(tags, diff.true, diff.col)
  }

  comment <- c(paste("# lib1:", lib1$libname, sep=" "),
               paste("# nseq1:", tagcount1, sep=" "),
               paste("# ntag1:", lib1$ntag, sep=" "),
               lib1$comment,
               paste("# lib2:", lib2$libname, sep=" "),
               paste("# nseq2:", tagcount2, sep=" "),
               paste("# ntag2:", lib2$ntag, sep=" "),
               lib2$comment)
  comp <- list(name=libname, ntag=length(tagnames), taglength=taglength, data=tags, comment=comment)
  class(comp) <- "sage.library.comparison"
  return(comp)
}

sage.test <- function(x, y, n1=sum(x), n2=sum(y))
#	Binomial probabilities for comparing SAGE libraries
#	Gordon Smyth
#	15 Nov 2003.  Last modified 21 Jan 2004.
{
	if(any(is.na(x)) || any(is.na(y))) stop("missing values not allowed")
	x <- round(x)
	y <- round(y)
	if(any(x<0) || any(y<0)) stop("x and y must be non-negative")
	if(length(x) != length(y)) stop("x and y must have same length")
	n1 <- round(n1)
	n2 <- round(n2)
	if(!missing(n1) && any(x>n1)) stop("x cannot be greater than n1")
	if(!missing(n2) && any(y>n2)) stop("y cannot be greater than n2")
	size <- x+y
	p.value <- rep(1,length(x))
	if(n1==n2) {
		i <- (size>0)
		if(any(i)) {
			x <- pmin(x[i],y[i])
			size <- size[i]
			p.value[i] <- pbinom(x,size=size,prob=0.5)+pbinom(size-x+0.5,size=size,prob=0.5,lower.tail=FALSE)
		}
		return(p.value)
	}
	prob <- n1/(n1+n2)
	if(any(big <- size>10000)) {
		ibig <- (1:length(x))[big]
		for (i in ibig) p.value[i] <- chisq.test(matrix(c(x[i],y[i],n1-x[i],n2-y[i]),2,2))$p.value
	}
	size0 <- size[size>0 & !big]
	if(length(size0)) for (isize in unique(size0)) {
		i <- (size==isize)
		p <- dbinom(0:isize,p=prob,size=isize)
		o <- order(p)
		cumsump <- cumsum(p[o])[order(o)]
		p.value[i] <- cumsump[x[i]+1]
	}
	p.value
}


# Graphics

difference.scatter.plot <- function(x1, y1, x2, y2, line.col="grey", col1="grey", col2="black", pch1=3, pch2=4, cex1=0.5, cex2=0.5,
                              xlim=range(x1,x2,na.rm=TRUE), ylim=range(range(y1,y2,na.rm=TRUE)), ...) {
  plot(0, 0, xlim=xlim, ylim=ylim, type="n", ...)
  lines(cbind(as.vector(rbind(x1,x2, NA)), as.vector(rbind(y1,y2, NA))), col=line.col)
  points(x1, y1, col=col1, pch=pch1, cex=cex1)
  points(x2, y2, col=col2, pch=pch2, cex=cex2)
}

plot.sage.library <- function(x, xlim=c(0,10), ...) {
  tt <- table(x$tags[,"count.raw"])
  tn <- as.integer(names(tt))
  yl <- max(tt, na.rm=TRUE)
  mmm<-max(tt)
  par(lab=c(10,5,7))
  plot(tn, tt, type="h", xlim=xlim, lty=1, lwd=2, xlab="tag abundance", ylab="frequency", ...)
  par(lab=c(5,5,7))
  if("count.adjusted" %in% names(x$tags)) {
    tt <- table(round(x$tags[,"count.adjusted"]))
    tn <- as.integer(names(tt))
    lines(tn+0.2, tt, type="h", col="gray20", lty=3, lwd=2)
  }
  if("count.true" %in% names(x$tags)) {
    tt <- table(x$tags[,"count.true"])
    tn <- as.integer(names(tt))
    lines(tn-0.2, tt, type="h", col="gray40", lty=2, lwd=2)
    legend(xlim[2], yl, c("observed counts n", "true counts m", "rounded estimated counts em"), lwd=2, lty=1:3, xjust=1, bty="n")
  } else {
    legend(xlim[2], yl, c("observed counts n", "rounded estimated counts em"), lwd=2, lty=c(1,3), xjust=1, bty="n")
  }
}

plot.sage.library.comparison <- function(x, xlab="A", ylab="M", main=x$name, pch=20, cex=0.5, ...) {
  A <- round(x$data$A, 2)
  M <- round(x$data$M.adjusted, 2)
  MA<-t(matrix(as.numeric(unlist(strsplit(unique(paste(A,M,sep=":")),split=":"))),nrow=2))
  plot(MA[,1], MA[,2], xlab=xlab, ylab=ylab, main=main, pch=pch, cex=cex, ...)
  abline(0,0)
}

# UTILITY functions

tagnum2tagmatrix <- function(tags, length) {
  t(matrix((rep(tags-1, each=length)%/% 4^((length-1):0)) %% 4, nrow=length))
}

tagmatrix2tagnum <- function(tags, length=ncol(tags)) {
  colSums(t(tags[,1:length])*4^((length-1):0))+1
}

tagnum2tagsequence <- function(tags, length) {
  new.tags <- t(matrix((rep(tags-1, each=length)%/% 4^((length-1):0)) %% 4, nrow=length))
  new.tags <- apply(new.tags, 1, paste, collapse="")
  chartr("0123", "acgt", new.tags)
}

tagsequence2tagnum <- function(tags, length) {
  new.tags <- tolower(unlist(strsplit(as.character(tags),"")))
  new.tags[!(new.tags=="a" | new.tags=="g" | new.tags=="c" | new.tags=="t" |
             new.tags=="s" | new.tags=="y" | new.tags=="b" | new.tags=="k")] <- "n"
  new.tags <- matrix(as.numeric(chartr("acgtnsybk", "012301112", new.tags)), nrow=length)
  colSums(new.tags*4^((length-1):0))+1
}

revcomp <- function(seq) {
  r <- lapply(strsplit(unlist(seq), ""), rev)
  r <- as.character(sapply(r, paste, collapse=""))
  chartr("ACGTUMRYKVHDBacgtumrykvhdb", "TGCAAKYRMBDHVtgcaakyrmbdhv", r)
}

summary.sage.library <- function(object, ...) {
  print.sage.library(object, ...)
  summary.factor(round(object$tags$count))
}

create.matrix.csr <- function(values, row.indices, col.indices, dim=NULL, eps = .Machine$double.eps) {
  o <- !(is.na(values) | is.na(row.indices) | is.na(col.indices) | as.numeric(values) < eps)
  ra <- as.numeric(values[o])
  ia <- as.integer(row.indices[o])
  ja <- as.integer(col.indices[o])
  if(is.null(dim)) dim <- c(max(ia), max(ja))
  library(SparseM)
  return(as.matrix.csr(new("matrix.coo", ra=ra, ia=ia, ja=ja, dimension=dim)))
}

table.sparse <- function(..., exclude = c(NA, NaN), dnn = list.names(...), deparse.level = 1) {
  # Modification of R function table. Omits all combinations that have count 0.
  # Instead of an array returns a list with vectors for dimnames and a vector for the counts.
  
  list.names <- function(...) {
    l <- as.list(substitute(list(...)))[-1]
    nm <- names(l)
    fixup <- if (is.null(nm)) 
      seq(along = l)
    else nm == ""
    dep <- sapply(l[fixup], function(x) switch(deparse.level + 
                                               1, "", if (is.symbol(x)) 
                                               as.character(x)
                                               else "", deparse(x)[1]))
    if (is.null(nm)) 
      dep
    else {
      nm[fixup] <- dep
      nm
    }
  }
  args <- list(...)
  if (length(args) == 0) 
    stop("nothing to tabulate")
  if (length(args) == 1 && is.list(args[[1]])) {
    args <- args[[1]]
    if (length(dnn) != length(args)) 
      dnn <- if (!is.null(argn <- names(args))) 
        argn
      else paste(dnn[1], 1:length(args), sep = ".")
  }
  bin <- 0
  lens <- NULL
  dims <- 1
  pd <- 1
  dn <- NULL
  res <- NULL
  for (a in args) {
    if (is.null(lens)) {
      lens <- length(a)
    } else if (length(a) != lens) {
      stop("all arguments must have the same length")
    }
    if (is.factor(a)) {
      cat <- a
    } else { cat <- factor(a, exclude = exclude) }
    nl <- length(l <- levels(cat))
    dims <- c(dims, nl)
    dn <- c(dn, list(l))
    bin <- bin + pd * (as.integer(cat) - 1)
    pd <- pd * nl
  }
  bin <- bin[!is.na(bin)]
  cat <- factor(bin)
  nl <- length(l <- levels(cat))
  bn <- as.integer(l)
  y <- tabulate(cat, nl)
  dims <- cumprod(dims)
  for (i in 1:length(args)) {
    res <- c(res, list((dn[[i]])[(bn %% dims[i+1] %/% dims[i])+1]))
  }
  res <- c(res, list(y))
  names(res) <- c(dnn, "count")
  class(res) <- "table.sparse"
  return(res)
}

print.sage.library <- function(x, ...) {
  write(c(paste("# libname:", x$libname, sep=" "),
          paste("# nseq:", x$nseq, sep=" "),
          paste("# ntag:", x$ntag, sep=" "),
          paste("# taglength:", x$taglength, sep=" "),
          x$comment, 
          "Fields:", paste(names(x), collapse=" "),
          "contents of field 'tags':", paste(names(x$tags), collapse=" "),
          "contents of field 'seqs':", paste(names(x$seqs), collapse=" ")
          ), "")
}

summary.sage.library.comparison <- function(object, ...) {
  print.sage.library.comparison(object, ...)
}

print.sage.library.comparison <- function(x, ...) {
  write(c(paste("# name:", x$name, sep=" "),
          paste("# ntag:", x$ntag, sep=" "),
          paste("# taglength:", x$taglength, sep=" "),
          x$comment, 
          "Fields:", paste(names(x), collapse=" "),
          "contents of field 'data':", paste(names(x$data), collapse=" "),
          ), "")
}
