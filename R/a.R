.First.lib <- function(libname,pkgname,where){
  if(.Platform$OS.type=="windows" && interactive() && .Platform$GUI=="Rgui"){
    addVigs2WinMenu("OrderedList")
  }  
}
 
