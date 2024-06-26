library("renv")
renv::init()
renv::snapshot()
renv::restore()

renv::deactivate()
renv::deactivate(clean = TRUE)
