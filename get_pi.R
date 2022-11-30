# update mixing weights pi based on logistic parameters
get_pi <- function(V,delta){
  exb = exp(V %*% delta)
  exb[which(is.infinite(exb))]=.Machine$double.xmax
  exb/apply(exb,1,sum)
}
