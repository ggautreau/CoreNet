setClass(Class="corenet",
         representation(dir_files="character", # Directory where the networks files per organism are located
                        union_graph="character", # name of an existing file 
                                                 # describing union of all networks
                        sep="character", # separator used in network files
                        nodes = "character", # vector containing a list of networks 
                        					 # nodes (reactions or genes)
                        organisms = "character", # vector containing 
                                                 # the name of all organisms                            
                        data="matrix",   # Binary data matrix each line is a node (reaction or gene), 
                                         # each column an organism. 1 is presence, 0 absence
                        neighbors="numeric", # Encoding of networks neighbors repeating 
                                         #('node number' 
                                         # 'number of neighbors' 
                                         # 'list of neighbors)
                        partitionList="list",
                        partitionOptions="list"                                       
                          ),
         prototype(dir_files=character(),
                   union_graph=character(),
                   sep=character(),
                   nodes=vector(mode="character",length=0),
                   organisms=vector(mode="character",length=0),
                   data=matrix(0),
                   neighbors=vector(mode="numeric",length=0),
                   partitionList=vector(mode="list"),
                   partitionOptions=vector(mode="list")))
         
#--------------------------------------------------------------
#  getAllNodes.corenet<- function(dir_files=".",union_graph='') 
#  return(list(AllNodes=AllNodes,edges=edges)) 
#  - AllNodes is a vector of character
#  - edges is a matrix with 2 columns describing the egdes  
#--------------------------------------------------------------  
       
getAllNodes.corenet<- function(dir_files=".",union_graph='',sep=" "){
if (union_graph==''){
# Find all the files related to all considered organism
 ListofFiles<-dir(dir_files)
 ListofSIFFiles<-NULL
 k<-0
 for (file in ListofFiles){
	strsplit(file,'\\.') -> splittedName
	if (length(splittedName[[1]])==2) {
		splittedName[[1]][2]->extension
	    splittedName[[1]][1]->baseName 
	    if (extension=="sif") {
		  k<-k+1;	ListofSIFFiles[k]<-file} 
   } 
 }
if (k==0){
	stop(paste("No .sif files in the directory",dir_files))
}

 # Creates the character vector of all nodes
 AllNodes <-  vector("character")
 for (file in ListofSIFFiles){
  fileID<-file(file.path(dir_files,file),blocking=TRUE)
	readLines(fileID,-1)->Nodes    
  close(fileID)
        AllNodes <- c(AllNodes,Nodes)
 }
AllNodes <- unique(AllNodes)
} else {  
          fileID<-file(file.path(dir_files,union_graph),blocking=TRUE)
          readLines(fileID,-1)->AllNodes    
          close(fileID)     
	}

strsplit(AllNodes,sep)->AllNodes

n<-length(AllNodes)
edges<-matrix(0,n,2)
k<-1
for (i in 1:n){
	if (length(AllNodes[[i]])>1){
		AllNodes[[i]] ->linkedNodes
    for (node in linkedNodes[-1]){
      edges[k,]<-sort(c(linkedNodes[1],node))      
		  k<-k+1
      }
		}
}
edges[-(k:n),]->edges    
unique(edges)->edges # get rid of the double

# All reactions, including isolated nodes
unique(unlist(AllNodes)) -> AllNodes
AllNodes[!(AllNodes==sep)] -> AllNodes
return(list(AllNodes=AllNodes,edges=edges))
}

#--------------------------------------------------------------
#  getNeighbors <- function(AllNodes,edges)
#  return(nei=unlist(neighbors)) 
#  Encoding  of neighbors repeating 
                                         #('node number' 
                                         # 'number of neighbors' 
                                         # 'list of neighbors)
#--------------------------------------------------------------  
getNeighbors.corenet <- function(AllNodes,edges){
# Convert edges numerical matrix (encoded using AllNodes)
edgesNum<-matrix(0,nrow(edges),2)
factor(AllNodes)->AllNodes
edgesNum[,1]<-as.numeric(factor(edges[,1],levels=AllNodes))
edgesNum[,2]<-as.numeric(factor(edges[,2],levels=AllNodes))

# find the neighbors of i
neighbors<-list()
for (i in 1:length(AllNodes)){
	neighborsOfi<-sort(c(edgesNum[edgesNum[,1]==i,2],edgesNum[edgesNum[,2]==i,1]))
	neighborsOfi<-list(c(i,length(neighborsOfi), neighborsOfi))
#	if (length(neighborsOfi[[1]])>2){
		neighbors<- append(neighbors,
		                   neighborsOfi)
#	} else {
    
#	}
    
}
return(nei=unlist(neighbors))  
}


#--------------------------------------------------------------
#  getPresence.corenet <-function(dir_files=".",union_graph='')
#  return(X) # X is a Binary data matrix 
#             # each line is a node (gene or metabolic reaction), 
#             # each column a species. 1 is presence, 0 absence
#--------------------------------------------------------------  
getPresence.corenet <-function(dir_files=".",union_graph='',sep=" "){
AllNodes <-  getAllNodes.corenet(dir_files,union_graph=union_graph,sep=sep)$AllNodes
ListofFiles<-dir(dir_files)
ListofSIFFiles<-NULL
k<-1
for (file in ListofFiles){
	if (!(union_graph=='')) 
	{ if (file != union_graph){
	strsplit(file,'\\.') -> splittedName
		if (length(splittedName[[1]])==2) {
	splittedName[[1]][2]->extension
	splittedName[[1]][1]->baseName 
	if (extension=="sif") {
		ListofSIFFiles[k]<-file; k<-k+1
	                       }
	                                       }
	                          }
	}
	else {
		strsplit(file,'\\.') -> splittedName
		if (length(splittedName[[1]])==2) {
			splittedName[[1]][2]->extension
	splittedName[[1]][1]->baseName 
	if (extension=="sif") {
		ListofSIFFiles[k]<-file; k<-k+1
	}
	}
	}
}
X<-matrix(0,length(AllNodes),length(ListofSIFFiles))
k<-1
for (file in ListofSIFFiles){

  fileID<-file(file.path(dir_files,file),blocking=TRUE)
  readLines(fileID,-1)->Nodes    
  close(fileID)
  
	unique(unlist(strsplit(Nodes,sep)))->Nodes
	# Let us get rid of separators
	Nodes[!(Nodes==sep)] -> Nodes
    X[,k]<-is.element(AllNodes,Nodes)
    k<-k+1
}
rownames(X) <- AllNodes
colnames(X) <- ListofSIFFiles
return(X)
}


#--------------------------------------------------------------
#  Methods for class corenet
#--------------------------------------------------------------  


setMethod("initialize",
	"corenet",
	function(.Object,dir_files=".",union_graph='',sep=" ") {
		  .Object@dir_files <- dir_files
		  .Object@union_graph <- union_graph
		  .Object@sep <- sep
		   nodes <- getAllNodes.corenet(dir_files,union_graph,sep=sep)
		  .Object@nodes <-    nodes$AllNodes
          .Object@data         <- getPresence.corenet(dir_files,union_graph,sep=sep)
          .Object@organisms <- colnames(.Object@data)
          .Object@neighbors <- getNeighbors.corenet(nodes$AllNodes,nodes$edges)
          return(.Object)
       })

setMethod("show",
	"corenet",
	 function(object) {
		cat("Object corenet builded from .sif files of dir",object@dir_files,"\n")
		if (object@union_graph=='') {cat("using no pre-given union graph \n")
			}	else {               cat("using union graph:",object@union_graph,"\n" )}	
	    numberOfSpecies<-ncol(object@data)
	    numberOfNodes<-nrow(object@data)
	    cat("Number of Species:",numberOfSpecies ,"\n")
	    cat("Number of Nodes:",numberOfNodes,"\n")		
	})

#--------------------------------------------------------------
## Neighborhing classification with EM algorithm
## Relies on an initial function written by P. Huppé for package Manor.
#--------------------------------------------------------------

nem.corenet <- function(X,nei,nk=NULL,nkMin=2,nkMax=2, family ="bernoulli", algo="nem",beta=1,bmod="optim",iters=2000,nbInit=30) {
    if (missing(X)) 
        stop("'X' is missing")
    if (missing(nei)) 
        stop("'nei' is missing")
    if (!any(bmod == c("fix", "optim"))) 
        stop("bmod must be fix or optim")
    if (!any(algo == c("ncem","nem","gem"))) 
        stop("Algo must be nem, ncem or gem") 
        
    familyNUM <- ( function(){
                         if (family=="bernoulli") return(1)
                         if (family=="normal") return(2) 
                        })()
    
     algoNUM <- ( function(){
                         if (algo=="nem") return(0)
                         if (algo=="ncem") return(1) 
                         if (algo=="gem") return(2)
                        })()
    
    
    if (bmod=="fix"){
	bmod.code<-1
	
        if(beta < 0)
            stop("beta must be >= 0")

        if(beta==0) {nei<-as.vector(rbind(seq(1:nrow(X)),rep(0,nrow(X)))); beta<-1}
     }	else {
     beta=1	
     bmod.code<-2	
     }
     if (!is.null(nk)){
    	nkMin=nk
    	nkMax=nk
     } 
    
    if ((nkMin<=0)|(nkMax<=0))
        stop("Number of classes must be > 0")
    if (nkMin>nkMax)
        stop("nkMin must be smaller than nkMax")

     
    if(iters < 100)
        stop("iters must be >= 100")

    NInds<-nrow(X)
    NVars<-ncol(X)
    lengthOfNei<-length(nei)
    span <-  nkMin:nkMax 
    NEMresult <- vector("list", length(span)) 
    #BICpM<-data.frame(nk=nkMin:nkMax,BICp=nkMin:nkMax) 
    q<-0
    for (nk in span){     
      q<-q+1  	
      res <- .C("nem",
              as.single(t(X)),
              as.integer(nei),
              as.integer(lengthOfNei),
              as.integer(NInds),
              as.integer(NVars), 
              as.integer(nk),
               as.integer(familyNUM),
              as.integer(algoNUM),
              as.single(beta),
              as.integer(bmod.code), 
              as.integer(iters),
              as.integer(nbInit),
              classification=integer(NInds),
              beta = single(1),
              mean = single(NVars*nk),
              Cik = single(NInds*nk),
              pseudoL = single(1),
              PACKAGE = "CoreExtension")
      pseudoL<-as.numeric(res$pseudoL)
      if (family =="bernoulli") {  dim2   <-NVars * nk + (nk-1)+(bmod.code -1)
                               } else {  dim2   <- (NVars+1 )*nk + (nk-1)+(bmod.code -1)}
   BICp <- 2*pseudoL - dim2*log(NInds) 
   NEMresult[[q]]$nk <- nk
   NEMresult[[q]]$BICp<-BICp
   print(q)
   NEMresult[[q]]$beta<- as.numeric(res$beta)
   NEMresult[[q]]$cluster<- res$classification
   NEMresult[[q]]$Cik<- matrix(res$Cik,NInds,nk,byrow=T)
  # if (nk==1) {
  # 	NEMresult[[q]]$mean<-colMeans(X)
  # }
  # else {
  	nbEmptyClasses<- nk-length(unique(res$classification))
    if  (nbEmptyClasses > 0) {
    	     EMPTYCLASS<-TRUE
    	     trueNk<-nk-nbEmptyClasses
    	     warning(paste(nbEmptyClasses, " empty class(es) when partitionning into ", nk, " classes."))
    	     } 
    	    else {
    	     EMPTYCLASS<-FALSE
    	     trueNk<-nk}
    if (trueNk==1) {effectifs<-sum(NEMresult[[q]]$Cik)   } else {effectifs<-colSums(NEMresult[[q]]$Cik)} 
    NEMresult[[q]]$mean<- t(NEMresult[[q]]$Cik / effectifs)%*%X
    #matrix(unlist(by(X,res$classification,colMeans)),trueNk,NVars,byrow=T)
  # }
      
   
   }
   
   
 return(list(output = NEMresult, 
             options=list(nk=nk,nkMin=nkMin,nkMax=nkMax,algo=algo,bmod=bmod,iters=iters))
           )
 }

setGeneric("nem",function(object,...) { standardGeneric("nem")})

setMethod("nem",
           "corenet",
           function(object,           	          
                      nk=NULL,
           	          nkMin=2,
                          nkMax=2,
                          family="bernoulli",
           	          algo="nem",
           	          beta=1,bmod=
           	          "optim",
           	          iters=2000,
           	          nbInit=5){         
           	result.nem<-nem.corenet(X=object@data,
           	          nei=object@neighbors,
           	          nk=nk,
           	          nkMin=nkMin,
           	          nkMax=nkMax,
           	          algo=algo,
           	          beta=beta,bmod=
           	          bmod,
           	          iters=iters,
           	          nbInit=nbInit) 	
           object@partitionList<-result.nem$output
           object@partitionOptions<-result.nem$options
           return(object)
           }
          
            )
    
setGeneric("getPartition",function(object,...) { standardGeneric("getPartition")})

setMethod("getPartition",
           "corenet",
           function(object,nk=2){
           	if (length(object@partitionList)==0) {
           		 stop("No partition is available: use nem(object) before extracting a partition...")}
           	partition<-NULL
           	 for (q in 1:length(object@partitionList)) {
           	   if (object@partitionList[[q]]$nk == nk) {       
                       partition<-object@partitionList[[q]]    
                       break;
                   } # end of  if (object@partitionList[[q]]$nk == nk)               
             } 	# end of for (q in 1:length(object@partitionList))
             if (is.null(partition)) {warning("No partition in the desired number of classes is available")}
             return(partition)
           }
          
            )

 
    
#-----------------------------------------------------
#
# plot 
#-----------------------------------------------------

    
            
setMethod("plot",
          "corenet",
          function(x,y,q=2,type="BICp",...){
          	  if ( length(x@partitionList) == 0) {
          	  	if (type=="mean"){
          	  	     image(t(x@data),col=gray(c(0,1)),axes = FALSE)
          	  	     }  else {
          	  		warning("No partition to compute BICp !!!, try plot(object,type='mean')")}
          	  	 	}
          	  else { 
          	  	if (type=="mean"){
          	  		n<-nrow(x@data)
          	  		numberOfPartitions<-length(x@partitionList)
          	  		for (j in 1:numberOfPartitions){
          	  		  dev.new()
          	  		  nk<-x@partitionList[[j]]$nk
          	  		  cluster<-x@partitionList[[j]]$cluster
         	  		  image(t(x@data[order(cluster),]),col=gray(c(0,1)),axes = FALSE,
         	  		        main=paste("Partition into",nk,"classes"))
         	  		  table(cluster)->limits # find the class limits
                     cumsum(limits)[1:(length(limits)-1)]+0.5->limits
                     abline(h=c(0,limits/n,1),col="red")

          	  		  }
          	  	}  else {
          	  	  #type =="BICp"
          	      x1<-unlist(lapply(x@partitionList,function(x) x$nk))
                  x2<-unlist(lapply(x@partitionList,function(x) x$BICp))
                  plot(x1,x2,xlab="Number of classes",ylab="BICp",main="Bayesian  Information Criterion computed with Pseudo-likelihood")
                 lines(x1,x2)
                 #abline(v=q,col="red",lty=2)
                   }
              }
          })            

		




         
         
         
