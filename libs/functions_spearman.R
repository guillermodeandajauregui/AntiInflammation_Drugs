#Spearman footrule

spearman_footrule <- function(ranks_i, ranks_j) {
                      distance = sum(abs(ranks_i - ranks_j))
                      return(distance)
}

#takes a ranked matrix, returns spearman footrule distance between columns
spearman_footrule_matrix <- function(matrix){
  return(apply(X = matrix, MARGIN = 2, FUN = function(i)
  {
    sapply(matrix, FUN = spearman_footrule, ranks_i = i)
  }
  )
  )
}