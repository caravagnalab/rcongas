.onLoad <- function(libname, pkgname)
{
  #  # =-=-=-=-=-=-
  #  # Required packages will be listed here
  #  # =-=-=-=-=-=-
  #  requirements = c(
  #    'pio',
  #    'easypar',
  #    'tidyverse',
  #    'tidygraph',
  #    'ggraph',
  #    'crayon',
  #    'igraph',
  #    'ggrepel',
  #    'RColorBrewer',
  #    'clisymbols',
  #    'entropy',
  #    'matrixcalc',
  #    'reshape2'
  #  )
  #
  #  suppressMessages(sapply(requirements, require, character.only = TRUE))

  # =-=-=-=-=-=-
  # Package options
  # =-=-=-=-=-=-
  options(pio.string_fg_colour = crayon::bgYellow$black)

  # =-=-=-=-=-=-
  # Header
  # =-=-=-=-=-=-

  rcongas_welcome_message =  getOption('rcongas_welcome_message', default = TRUE)

  if (rcongas_welcome_message)
  {
    pk = 'Rcongas'
    pk_l = 'Copy-Number genotyping from single cells'
    www = "https://militeee.github.io/Rcongas/"
    em = "xxxx@gmail.com"

    cli::cli_alert_success(
      'Loading {.field {pk}}, {.emph \'{pk_l}\'}. Support : {.url { www}}' )

    options(rcongas_welcome_message = FALSE)
  }

  invisible()
}
