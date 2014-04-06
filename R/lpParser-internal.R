lpParser <- setClass(
  'lpParser', 
  slots = c(indexDim = 'numeric', 
    indexLength = 'numeric', 
    indexNames = 'character', 
    indexList = 'list',
    constDim = 'numeric',
    constNames = 'character',
    constList = 'list',
    constIndex = 'list',
    varDim = 'numeric',
    varLength = 'numeric',
    varNames = 'character',
    varIndex = 'list',
    varList = 'list',
    eqDim = 'numeric',
    eqLength = 'numeric',
    eqNames = 'character',
    # eqNamesIndex = 'character',
    eqList = 'list',
    eqIndex = 'list',
    eqRhs = 'numeric',
    eqLhs = 'matrix',
    eqType = 'numeric',
    obj = 'dsparseVector',
    lpModel = 'lpExtPtr',
    status = 'numeric',
    objective = 'numeric',
    varVector = 'numeric',
    varArrays = 'list',
    boundsDim = 'numeric',
    boundsList ='list',
    sense = 'character',
    type = 'numeric'
  ),
  prototype = list(indexDim = 1, 
   indexLength = 1, 
   indexNames = 'single', 
   indexList = list(''),
   constDim = 0,
   constNames = NULL,
   constList = list(),
   constIndex= list(),
   varDim = 0,
   varLength = 0,
   varNames = NULL,
   varIndex = list(),
   varList = list(),
   eqDim = 0,
   eqLength = 0,
   eqNames = NULL,
   # eqNamesIndex = NULL,
   eqIndex = list(),
   eqList = list(),
   eqRhs = NULL,
   eqLhs = NULL,
   eqType = NULL,
   obj = NULL,
   lpModel = NULL,
   status = NULL,
   objective = NULL,
   varVector = NULL,
   varArrays = NULL,
   boundsDim = 0,
   boundsList = NULL,
   sense = 'min',
   type = c('<=' = 1, '=<' = 1, '<' =1,
    '>=' = 2, '=>' =  2, '>' = 2,
    '=' = 3, '==' = 3))
)


.stopifnot <- function (..., text) {
  n <- length(ll <- list(...))
  if (n == 0L) 
    return(invisible())
  mc <- match.call()
  for (i in 1L:n) if (!(is.logical(r <- ll[[i]]) && !any(is.na(r)) && 
                          all(r))) {
    ch <- deparse(mc[[i + 1]], width.cutoff = 60L)
    if (length(ch) > 1L) 
      ch <- paste(ch[1L], "....")
    stop(paste(text,sprintf(ngettext(length(r), "%s is not TRUE", "%s are not all TRUE"), 
                            ch), sep='\n'), call. = FALSE, domain = NA)
  }
  invisible()
}

.message = c("optimal solution found",
             "the model is sub-optimal",
             "the model is infeasible",
             "the model is unbounded",
             "the model is degenerate",
             "numerical failure encountered",
             "process aborted",
             "timeout",
             "the model was solved by presolve",
             "the branch and bound routine failed",
             "the branch and bound was stopped because of a break-at-first or break-at-value",
             "a feasible branch and bound solution was found",
             "no feasible branch and bound solution was found")

.stSplit = function(string){
  sign = 1
  if (length(.minus(string)))  sign = -1
  string = .rmPM(string)
  if (0 < length(grep("[*]", string))) {
    a = strsplit(string, '[\\*]')[[1]]
    fac = gsub('[ ]','',a[1])
    var = gsub('[ ]','',a[2])
  } else {
    fac = '1'
    var = gsub('[ ]','',string)
  }
  return(c(sign = sign, factor = fac, variable = var))
}


.addSparseVectors <- function(x){
  y = x[[1]]
  n = length(x)
  if(n > 1) for (i in seq(2,n)) y = y + x[[i]]
  return(y)
}

.replace <- function(x, oldSlice, newSlice){
  sapply(x, function(x) newSlice[match(x,oldSlice)])
}

.arrayReplace = function(array,vector){
  D <- dim(array)
  y <- vector[array]
  dim(y) <- D
  return(y)
}

.checkSolution = function(object){
  if (is.null(object@status)) {
    print('The model has to be solved first!')
    return(FALSE)
  }
  if (object@status != 0) {
    print('The model has not been solved sucsefully')
    return(FALSE)
  }
  return(TRUE)
}

setGeneric('.makeMatrixPos',function(object, ...){standardGeneric('.makeMatrixPos')})

setMethod(f = '.makeMatrixPos',
  signature = 'lpParser',
  definition = function(object, grid, sliceList, factItem, env){
    with(as.list(factItem), {
      y = object@eqLength + 1:nrow(grid)
      # Constant
      constNr = match(factor, object@constNames)
      if (!is.na(constNr)){
        # gridPos gives the position in the grid of intersection between eqIndex og constIndex
        gridPos = which(object@eqIndex[[object@eqDim]] %in% object@constIndex[[constNr]])
        # constPos gives the position in the constIndex which intersec with eqIndex and therfore is in grid
        constPos = which(object@constIndex[[constNr]] %in% object@eqIndex[[object@eqDim]][gridPos])
        # arrayPos gives a list of seqenses coresonding to contsIndex
        arrayPos = lapply(object@indexLength[object@constIndex[[constNr]]], seq)
        # For each line in the grid the values of the grid is replaced with the seq in the arrayPos
        # and the proper values retrived from the list
        val = apply(grid[,gridPos,drop = FALSE], 1, function(x) {
          arrayPos[constPos] <- x
          as.vector(do.call('[', c(list(object@constList[[constNr]]), arrayPos)))
        }) * as.numeric(sign)
        if (length(dim(val)) == 0){
          dim(val) <- c(1,length(val))
        }
      } else {
        .stopifnot(length(grep('[0-9]',factor)) & length((grep('[^0-9]',factor)))==0, 
                   text = paste('Left hand side factor:', factor, '- not recornised 2'))
        val = rep(eval(parse(text = factor)), length.out = nrow(grid)) * as.numeric(sign)
        dim(val) <- c(1,length(val))
      }
      # variable
      varStr = .getString(variable)
      varSlice = .getSlice(variable)
      varNr = match(varStr, object@varNames)
      .stopifnot(!is.na(varNr), text = paste('Left hand variable:', variable, 'not recornsied'))
      if (all(object@varIndex[[varNr]] == 1)){ #variable is 'single'
        xPos = rep(object@varList[[varNr]], length.out = nrow(grid))
        dim(xPos) = c(1,length(xPos))
      } else {
        # gridPos gives the position in the grid of intersection between eqIndex og constIndex
        gridPos = which(object@eqIndex[[object@eqDim]] %in% object@varIndex[[varNr]])
        if (!length(gridPos)){
          # CHECK FOR INDEX
          # if the grid pos is empty, ei no intersection between eqIndex og constIndex,
          # the variable possition is return as a whole
          xPos = apply(grid[,,drop = FALSE], 1, function(x) as.vector(object@varList[[varNr]]))
        } else {
          # modify grid by replcing index with the given slice
          # CHECK LENGTH OF INDEX
          varGrid <- grid
          indexSlicePos = which(!is.na(sliceList))
          if (!is.na(varSlice)) varGrid[,indexSlicePos] <- .replace(
            varGrid[,indexSlicePos], 
            eval(parse(text = sliceList[[indexSlicePos]]), envir = env),
            eval(parse(text = varSlice), envir = env))
          # varPos gives the position in the varIndex which intersec with eqIndex and therfore is in grid
          varPos = which(object@varIndex[[varNr]] %in% object@eqIndex[[object@eqDim]][gridPos])
          
          # arrayPos gives a list of seqenses coresonding to varIndex
          arrayPos = lapply(object@indexLength[object@varIndex[[varNr]]], seq)
          # For each line in the grid the values of the grid is replaced with the seq in the arrayPos
          # and the proper values retrived from the list
          xPos = apply(varGrid[,gridPos,drop = FALSE], 1, function(x) {
            arrayPos[varPos] <- x
            as.vector(do.call('[',c(list(object@varList[[varNr]]),arrayPos)))
          })
          if (is.null(dim(xPos))) dim(xPos) = c(1,length(xPos))
        }
      }
      n = nrow(xPos)
      .stopifnot(nrow(val)<=n, 
                 text = paste('Here is some problem: A constant:', factor, 'with dim', nrow(val),
                              'have bigger dimensions than the variable:', variable, 'with dim', n))
      cbind(x = as.vector(xPos),y = rep(y, each=n),val = as.vector(apply(val,2,rep, length.out=n)))
    })
  }
)




setGeneric('.getConstVal',function(object, ...){standardGeneric('.getConstVal')})

setMethod(f = '.getConstVal',
          signature = 'lpParser',
          definition = function(object, grid, constNr){
            if (length(object@constList[[constNr]]) == 1) {
              return(rep(object@constList[[constNr]], length.out = nrow(grid)))
            }
            .stopifnot(object@constIndex[[constNr]] %in% object@eqIndex[[object@eqDim]],
                       text = 'Constant index do not match equation index')
            gridPos = match(object@constIndex[[constNr]], object@eqIndex[[object@eqDim]])
            apply(grid[,gridPos,drop=F],1,function(x) do.call('[',c(list(object@constList[[constNr]]), x)))
          }
)

setGeneric('.makeObjItems',function(object, ...){standardGeneric('.makeObjItems')})

setMethod(f = '.makeObjItems',
          signature = 'lpParser',
          definition = function(object, factItem){
            with(as.list(factItem), {
              objItem = sparseVector(x=0, i =1, length=object@varLength)
              # variable
              varNr = match(variable, object@varNames)
              .stopifnot(!is.na(varNr), text = paste('Left hand variable:', variable, 'not recornsied'))
              a = cbind(index = object@varIndex[[varNr]], var = dim(object@varList[[varNr]]), const = 1)
              if (nrow(a)==1) a = rbind(a, c(1,1,1))
              # Constant
              constNr = match(factor, object@constNames)
              if (!is.na(constNr)){
                # objPos gives the position in the objVarList of intersection between objIndex og constIndex
                objPos = which(object@varIndex[[varNr]] %in% object@constIndex[[constNr]])
                
                constList = object@constList[[constNr]] * as.numeric(sign)
                a[objPos,3] = a[objPos,2]
              } else {
                .stopifnot(length(grep('[0-9]',factor)) & length((grep('[^0-9]',factor)))==0, 
                           text = paste('Constant factor:', factor, '- not recornised 2'))
                constList = eval(parse(text = factor)) * as.numeric(sign)
              }
              dim(constList) = a[,3]
              objItem[as.vector(object@varList[[varNr]])] = as.vector(.expandArray(constList, a[,2]))
              return(objItem)
            })
          }
)

.expandArray <- function(x, newDim){
  oldDim <- dim(x)
  slice = mapply(function(old,new){
    if (old==new) return(seq(old))
    if (old==1) return(rep(1,new))
    .stopifnot(FALSE, 'Dimension mismatch in .expandArray')},
    oldDim, newDim)
  do.call("[",c(list(x),slice))
}

setGeneric('.makeObjItems',function(object, ...){standardGeneric('.makeObjItems')})

setMethod(f = '.makeObjItems',
          signature = 'lpParser',
          definition = function(object, factItem){
            with(as.list(factItem), {
              objItem = sparseVector(x=0, i =1, length=object@varLength)
              # variable
              varNr = match(variable, object@varNames)
              .stopifnot(!is.na(varNr), text = paste('Left hand variable:', variable, 'not recornsied'))
              a = cbind(index = object@varIndex[[varNr]], var = dim(object@varList[[varNr]]), const = 1)
              if (nrow(a)==1) a = rbind(a, c(1,1,1))
              # Constant
              constNr = match(factor, object@constNames)
              if (!is.na(constNr)){
                # objPos gives the position in the objVarList of intersection between objIndex og constIndex
                objPos = which(object@varIndex[[varNr]] %in% object@constIndex[[constNr]])
                
                constList = object@constList[[constNr]] * as.numeric(sign)
                a[objPos,3] = a[objPos,2]
              } else {
                .stopifnot(length(grep('[0-9]',factor)) & length((grep('[^0-9]',factor)))==0, 
                           text = paste('Constant factor:', factor, '- not recornised 2'))
                constList = eval(parse(text = factor)) * as.numeric(sign)
              }
              dim(constList) = a[,3]
              objItem[as.vector(object@varList[[varNr]])] = as.vector(.expandArray(constList, a[,2]))
              return(objItem)
            })
          }
)

.expandArray <- function(x, newDim){
  oldDim <- dim(x)
  slice = mapply(function(old,new){
    if (old==new) return(seq(old))
    if (old==1) return(rep(1,new))
    .stopifnot(FALSE, 'Dimension mismatch in .expandArray')},
    oldDim, newDim)
  do.call("[",c(list(x),slice))
}

.getParenArg <- function(string){
  sSplit <- strsplit(string[1],'')[[1]]
  substr(string, min(grep('[(]',sSplit)) + 1,
         max(grep('[)]',sSplit))-1)
}

.getString = function(str) {
  strsplit(str,'[[]')[[1]][1]
}

.getSlice = function(str){
  strsplit(strsplit(str,'[[]')[[1]][2],'[]]')[[1]][1]
}

.pmPos <- function(string){
  sSplit = strsplit(string,'')[[1]]
  setdiff(grep('[-+]',strsplit(string,'')[[1]]),
          unlist(apply(cbind(grep('[[]',sSplit),grep('[]]',sSplit)),
                       1, function(x) x[1]:x[2])))
}

.rmPM <- function(string){
  substr(string,max(0,.pmPos(string)+1),1000L)
}

.minus <- function(string){
  sSplit = strsplit(string,'')[[1]]
  setdiff(grep('[-]',strsplit(string,'')[[1]]),
          unlist(apply(cbind(grep('[[]',sSplit),grep('[]]',sSplit)),
                       1, function(x) x[1]:x[2])))
}