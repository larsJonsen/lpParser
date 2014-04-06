
setGeneric("Equation<-",function(object, ...){standardGeneric("Equation<-")})

setReplaceMethod(
  f="Equation",
  signature="lpParser",
  definition=function(object, ...){
    eqList <- (...)
    env = parent.frame()
    .stopifnot(is.character(eqList), text = 'Eqations have to be characters')
    eqSplit <- strsplit(eqList,';')
    .stopifnot(sapply(eqSplit,function(x) length(x)==2),
              text = 'Each equation need one only semi colon')
    #env <- environment()
    for (eqString in eqSplit){
      object@eqDim = object@eqDim +1
      # Index
      indexRaw <- unlist(
        lapply(strsplit(.getParenArg(eqString[1]), ','),
          function(x) gsub(' ','',x)))
      indexList <- lapply(indexRaw, .getString)
      sliceList <- lapply(indexRaw, .getSlice)
      .stopifnot(indexList %in% object@indexNames, 
                text = 'Indexs in equation are not declared as indexs')
      
      object@eqIndex[[object@eqDim]] = sort(match(indexList, object@indexNames))
      ## Index grid
      # Index grid gives all combinationions of the indexes in a grid
      grid = expand.grid(mapply(function(x,y){
          if (is.na(x)) return(seq(y))
          return(eval(parse(text = x), envir = env))
        },
        sliceList,
        object@indexLength[object@eqIndex[[object@eqDim]]]))
      # Numbers of equations ei length of grid
      eqNr = object@eqLength + 1:nrow(grid)
            # equation name
      eqNam <- strsplit(eqString[1],'[(]')[[1]][1]
      object@eqNames = c(object@eqNames, eqNam)
      namGrid = apply(grid,1,
          function(x) mapply(function(x,y) y[x],
            as.list(x), object@indexList[object@eqIndex[[object@eqDim]]]))
      if (length(dim(namGrid))==0) dim(namGrid) = c(1,length(namGrid))
      namGrid = t(namGrid)
#       nameList = apply(namGrid, 1, 
#                        function (x) paste(eqNam, paste(Map(as.character,x),collapse = '_'),sep ='_'))
#       object@eqNamesIndex = c(object@eqNamesIndex, gsub(' ','_',nameList))
      # equation index grid
      namGrid = as.data.frame(namGrid)
      names(namGrid) <- object@indexNames[object@eqIndex[[object@eqDim]]]
      rownames(namGrid) <- eqNr
      object@eqList[[eqNam]] = namGrid
      # equation type
      a <- grep('[<=>]', strsplit(eqString[2],'')[[1]])
      typeStr <- substr(eqString[2], min(a), max(a))
      .stopifnot(typeStr %in% names(object@type), 
                text ='Eqution sign not recornised \n use <=. >= or ==')
      object@eqType = c(object@eqType,  rep(object@type[typeStr], length.out = nrow(grid)))
      # right hand side
      rhsStr = gsub(' ', '', substr(eqString[2], max(a) + 1, 1000L))
      ## CHECK FOR [] IN RHS
      constNr = match(rhsStr ,object@constNames)
      if (!is.na(constNr)){
        object@eqRhs = c(object@eqRhs, .getConstVal(object, grid, constNr))
      } else {
        .stopifnot(is.numeric(eval(parse(text = rhsStr))), 
                  text = 'Right hand side not constant or numeric')
        object@eqRhs = c(object@eqRhs, rep(eval(parse(text = rhsStr)), length.out = nrow(grid))) 
      }
      # left hand side
      lhsStr = substr(eqString[2], 0, min(a)-1)
      p = .pmPos(lhsStr)
      pp = cbind(c(0,p),c(p-1,1000L))
      termList = apply(pp,1,function(x) substr(lhsStr,x[1],x[2]))
      
      factorList <- lapply(termList, .stSplit)
      MatrixPos = lapply(factorList,function(x) .makeMatrixPos(object, grid, sliceList, x, env))
      object@eqLhs = do.call(rbind, c(list(object@eqLhs), MatrixPos))
      # update total number of eqations
      object@eqLength = object@eqLength + nrow(grid) 
    } 
    return (object)
  }
)


setGeneric("Variable<-",function(object, ...){standardGeneric("Variable<-")})

setReplaceMethod(
  f="Variable",
  signature="lpParser",
  definition=function(object, ...){
    var <- (...)
    Vname <- names(var)
    n <- length(var)
    for (i in seq(n)){
      .stopifnot(var[[i]] %in% object@indexNames, 
                text = 'Variable index not in the declared indexs') 
      object@varDim = object@varDim + 1
      object@varNames = c(object@varNames, Vname[i])
      object@varIndex[[object@varDim]] = sort(match(var[[i]], object@indexNames))
      object@varList[[object@varDim]] = object@varLength +
        array(1:prod(object@indexLength[object@varIndex[[object@varDim]]]),
              object@indexLength[object@varIndex[[object@varDim]]])
      
      namGrid = expand.grid(object@indexList[object@varIndex[[object@varDim]]])
      nameList = apply(namGrid,1, function(x) paste(Vname[i],paste(x,collapse = '_'),sep='_'))
#       object@varNamesIndex = c(object@varNamesIndex,gsub(' ','_',nameList))
      
      object@varLength = object@varLength + 
        prod(object@indexLength[object@varIndex[[object@varDim]]])
    }
    return (object)
  }
)


setGeneric("Index<-",function(object, ...){standardGeneric("Index<-")})
 
setReplaceMethod(
  f="Index",
  signature="lpParser",
  definition=function(object, ...){
    index <- (...)
    Iname <- names(index)
    n <- length(index)
    for (i in seq(n)){
      object@indexDim = object@indexDim +1
      object@indexLength = c(object@indexLength, length(index[[i]]))
      object@indexNames = c(object@indexNames, Iname[i])
      object@indexList[[object@indexDim]] = index[[i]]}
    return (object)
    }
)

setGeneric("Constant<-",function(object, ...){standardGeneric("Constant<-")})

setReplaceMethod(
  f="Constant",
  signature="lpParser",
  definition=function(object, ...){
    const <- (...)
    Cname <- names(const)
    n <- length(const)
    for (i in seq(n)){
      .stopifnot(const[[i]] %in% object@indexNames,
                text = 'Constant index not in declared indexs ') 
      object@constDim = object@constDim + 1
      object@constNames = c(object@constNames, Cname[i])
      object@constIndex[[object@constDim]] = sort(match(const[[i]], object@indexNames))
      object@constList[[object@constDim]] = array(NA, object@indexLength[object@constIndex[[object@constDim]]])
    }
    return (object)
  }
)

setGeneric("setConstant<-",function(object, ...){standardGeneric("setConstant<-")})

setReplaceMethod(
  f="setConstant",
  signature="lpParser",
  definition=function(object, ...){
    const <- (...)
    Cname <- names(const)
    n <- length(const)
    for (i in seq(n)){
      .stopifnot(Cname[i] %in% object@constNames, 
                text = 'Constant name not in the declared cosntants')
      j <- match(Cname[i], object@constNames)
      .stopifnot(dim(const[[i]]) == dim(object@constList[[j]]),
                text = 'Constant diensions do not match the declared indexs dim') 
      object@constList[[j]] <- const[[i]]
    }
    return (object)
  }
)


  
setGeneric('lhsMatrix',function(object, ...){standardGeneric('lhsMatrix')})

setMethod(f = 'lhsMatrix',
          signature = 'lpParser',
          definition = function(object){
            require(Matrix)
            return(sparseMatrix(j = object@eqLhs[,1], i = object@eqLhs[,2], x= object@eqLhs[,3], 
                                dims = c(object@eqLength, object@varLength)))
          }
)

setGeneric('lhsImage',function(object, ...){standardGeneric('lhsImage')})

setMethod(f = 'lhsImage',
          signature = 'lpParser',
          definition = function(object){
            require(Matrix)
            image(lhsMatrix(object))
          }
)

setGeneric("Object<-",function(object, ...){standardGeneric("Object<-")})

setReplaceMethod(
  f="Object",
  signature="lpParser",
  definition=function(object, ...){
    objStr <- (...)
    .stopifnot(is.character(objStr), text = 'Objects have to be characters')
    p = grep('[-+]',strsplit(objStr,'')[[1]])
    pp = cbind(c(0,p),c(p-1,1000L))
    termList = apply(pp,1,function(x) substr(objStr,x[1],x[2]))
    factorList <- lapply(termList, .stSplit)
    objList = lapply(factorList,function(x) .makeObjItems(object, x))
    object@obj = .addSparseVectors(objList)
    return(object)
  }
)

setGeneric("Bounds<-",function(object, ...){standardGeneric("Bounds<-")})

setReplaceMethod(
  f="Bounds",
  signature="lpParser",
  definition=function(object, ...){
    argList <- (...)
    argNames = names(argList)
    with(argList, {
      if (!('variable' %in% argNames)) variable = NULL
      if (!('lower' %in% argNames)) lower = NULL
      if (!('upper' %in% argNames)) upper = NULL
      if (!('columns' %in% argNames)) columns = NULL
    
      if (!is.null(variable)){
        .stopifnot(is.null(columns), text = ' variables and colunmns can not be set at same time')
        varNr = match(variable, object@varNames)
        .stopifnot(!is.na(varNr), text = 'Variable not recornised')
        columns = as.vector(object@varList[[1]])
        
      } else {
        if (is.null(columns)) columns = 1:object@varLength
      }
      object@boundsDim = object@boundsDim + 1
      object@boundsList[[object@boundsDim]] = list(
        lower = suppressWarnings(rep(lower, length.out = length(columns))),
        upper = suppressWarnings(rep(upper, length.out = length(columns))),
        columns = columns)
      return(object)
    })
  }
)

setGeneric("Solve",function(object, ...){standardGeneric("Solve")})

setMethod(
  f="Solve",
  signature="lpParser",
  definition=function(object, ...){
    env = parent.frame()
    objectName <- as.character(as.list(match.call())[['object']]) 
    object@lpModel = make.lp(nrow = object@eqLength, ncol = object@varLength)
    set.rhs(object@lpModel, object@eqRhs)
    set.constr.type(object@lpModel, object@eqType)
    con = lhsMatrix(object)
    for (i in seq(object@varLength)){
      set.column(object@lpModel, i,
                 con@x[(con@p[i]+1):con@p[i+1]],
                 indices=con@i[(con@p[i]+1):con@p[i+1]]+1)
    }
    set.objfn(object@lpModel, object@obj@x, object@obj@i)
    do.call(set.bounds, c(list(object@lpModel),Bounds(object)))
    lp.control(object@lpModel, sense= object@sense)
    object@status = solve(object@lpModel)
    cat(paste( 'Status rescived from lpSolve:', object@status, '\n ',
               .message[object@status + 1]))
    if (object@status==0){
      primal = get.primal.solution(object@lpModel)
      object@objective = primal[1]
      object@varVector = primal[(2+object@eqLength):length(primal)]
      eqVal = primal[2:(1+object@eqLength)]
      object@varArrays = lapply(object@varList, 
                                .arrayReplace, object@varVector)
      names(object@varArrays) = object@varNames
      for (i in seq(along.with = object@varIndex)){
        namList = lapply(object@varIndex[[i]], function(x){
          return(object@indexList[[x]])})
        names(namList) = object@indexNames[object@varIndex[[1]]]
        dimnames(object@varArrays[[i]]) = namList}
      
      dual = get.dual.solution(object@lpModel)
      eqLambda = dual[2: (1+object@eqLength)]
      for (i in seq(along = object@eqList)){
        object@eqList[[i]]$value = eqVal[as.numeric(rownames(object@eqList[[i]]))]
        object@eqList[[i]]$lambda = eqLambda[as.numeric(rownames(object@eqList[[i]]))]
      }
    }
    # modify the origonal object 
    # NOT THE R WAY !!
    assign(objectName, object, envir = env)
    # return the object invisible  
    # for security if the output i assigned the R way
    invisible(object)
  }
)


setGeneric("getVariables",function(object, ...){standardGeneric("getVariables")})

setMethod(
  f="getVariables",
  signature="lpParser",
  definition=function(object, ...){
    if (.checkSolution(object)) return(object@varArrays)
    invisible(NULL)
  }
)

setGeneric("getConstraints",function(object, ...){standardGeneric("getConstraints")})

setMethod(
  f="getConstraints",
  signature="lpParser",
  definition=function(object, ...){
    if (.checkSolution(object)) return(object@eqList)
    invisible(NULL)
  }
)

setGeneric("getObjective",function(object, ...){standardGeneric("getObjective")})

setMethod(
  f="getObjective",
  signature="lpParser",
  definition=function(object, ...){
    if (.checkSolution(object)) return(get.objective(object@lpModel))
    invisible(NULL)
  }
)
setGeneric("Bounds",function(object, ...){standardGeneric("Bounds")})

setMethod(
  f="Bounds",
  signature="lpParser",
  definition=function(object, ...){
    if(!is.null(object@status)) if(object@status == 0) return(get.bounds(object@lpModel))
    lower = rep(0, object@varLength)
    upper = rep(Inf, object@varLength)
  lapply(object@boundsList, function(x){
    if (!is.null(x$lower)) lower[x$columns] <<- x$lower
    if (!is.null(x$upper)) upper[x$columns] <<- x$upper
  })
  return(list(lower = lower, upper = upper))
  }
)

setGeneric("Sense",function(object, ...){standardGeneric("Sense")})

setMethod(
  f="Sense",
  signature="lpParser",
  definition=function(object, ...){
    env = parent.frame()
    objectName <- as.character(as.list(match.call())[['object']])    
    if (missing(...)) {
      print(object@sense)
    } else {
      x <- (...)
      .stopifnot(x %in% c('min','max'),
              text ="Sense expect either 'min' or 'max' ")
      object@sense = x
       # modify the origonal object 
      # NOT THE R WAY !!
      assign(objectName, object, envir = env)
    }
    # return the object invisible  
    # for security if the output i assigned the R way
    invisible(object)
  }
)