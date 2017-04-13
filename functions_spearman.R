#Spearman footrule

spearman_footrule <- function(ranks_i, ranks_j) {
                      distance = sum(abs(ranks_i - ranks_j))
                      return(distance)
                      }