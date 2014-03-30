\name{lpParser-class}
\Rdversion{1.1}
\docType{class}
\alias{lpParser-class}
%\alias{lpParser}
%\alias{.getConstVal,lpParser-method}
%\alias{.makeMatrixPos,lpParser-method}
%\alias{.makeObjItems,lpParser-method}
%\alias{Bounds,lpParser-method}
%\alias{Bounds<-,lpParser-method}
%\alias{Constant<-,lpParser-method}
%\alias{Equation<-,lpParser-method}
%\alias{getConstraints,lpParser-method}
%\alias{getVariables,lpParser-method}
%\alias{Index<-,lpParser-method}
%\alias{lhsImage,lpParser-method}
%\alias{lhsMatrix,lpParser-method}
%\alias{Object<-,lpParser-method}
%\alias{Sense,lpParser-method}
%\alias{setConstant<-,lpParser-method}
%\alias{Solve,lpParser-method}
%\alias{Variable<-,lpParser-method}

\title{Class \code{"lpParser"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{lpParser()}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{indexDim}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{indexLength}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{indexNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{indexList}:}{Object of class \code{"list"} ~~ }
    \item{\code{constDim}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{constNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{constList}:}{Object of class \code{"list"} ~~ }
    \item{\code{constIndex}:}{Object of class \code{"list"} ~~ }
    \item{\code{varDim}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{varLength}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{varNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{varIndex}:}{Object of class \code{"list"} ~~ }
    \item{\code{varList}:}{Object of class \code{"list"} ~~ }
    \item{\code{eqDim}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{eqLength}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{eqNames}:}{Object of class \code{"character"} ~~ }
    \item{\code{eqList}:}{Object of class \code{"list"} ~~ }
    \item{\code{eqIndex}:}{Object of class \code{"list"} ~~ }
    \item{\code{eqRhs}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{eqLhs}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{eqType}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{obj}:}{Object of class \code{"dsparseVector"} ~~ }
    \item{\code{lplpParser}:}{Object of class \code{"lpExtPtr"} ~~ }
    \item{\code{status}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{objective}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{varVector}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{varArrays}:}{Object of class \code{"list"} ~~ }
    \item{\code{boundsDim}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{boundsList}:}{Object of class \code{"list"} ~~ }
    \item{\code{sense}:}{Object of class \code{"character"} ~~ }
    \item{\code{type}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{.getConstVal}{\code{signature(object = "lpParser")}: ... }
    \item{.makeMatrixPos}{\code{signature(object = "lpParser")}: ... }
    \item{.makeObjItems}{\code{signature(object = "lpParser")}: ... }
    \item{Bounds}{\code{signature(object = "lpParser")}: ... }
    \item{Bounds<-}{\code{signature(object = "lpParser")}: ... }
    \item{Constant<-}{\code{signature(object = "lpParser")}: ... }
    \item{Equation<-}{\code{signature(object = "lpParser")}: ... }
    \item{getConstraints}{\code{signature(object = "lpParser")}: ... }
    \item{getVariables}{\code{signature(object = "lpParser")}: ... }
    \item{Index<-}{\code{signature(object = "lpParser")}: ... }
    \item{lhsImage}{\code{signature(object = "lpParser")}: ... }
    \item{lhsMatrix}{\code{signature(object = "lpParser")}: ... }
    \item{Object<-}{\code{signature(object = "lpParser")}: ... }
    \item{Sense}{\code{signature(object = "lpParser")}: ... }
    \item{setConstant<-}{\code{signature(object = "lpParser")}: ... }
    \item{Solve}{\code{signature(object = "lpParser")}: ... }
    \item{Variable<-}{\code{signature(object = "lpParser")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("lpParser")
}
\keyword{classes}
\usage{
lpParser(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
