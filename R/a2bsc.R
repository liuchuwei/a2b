require(pryr)

# DimRe class
setClass("DimRe", slots = list(umap = "data.frame",
                               tnse = "data.frame"))

# Annotation class
setClass("Ann", slots = list(ann = "data.frame",
                             markers = "data.frame"))


setClass("Annotation", slots = list(main = "Ann",
                             Tcell = "Ann",
                             Bcell = "Ann",
                             DCcell = "Ann",
                             MACcell = "Ann",
                             NKcell = "Ann"))


# a2bsc class
setClass("a2bsc",slots=list(Exp="Matrix",
                            pData="data.frame",
                            fData = "data.frame",
                            dimRe = "DimRe",
                            Ann = "Annotation"))
