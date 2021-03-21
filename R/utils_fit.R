tensorize <-  function(x){

  torch <- reticulate::import("torch")
  ret <-  torch$tensor(x)
  return(ret)
}


assemble_input <- function(x){

}
