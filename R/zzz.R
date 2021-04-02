.onLoad <- function(libname, pkgname)
{

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
    pk = '(R)CONGAS+'
    pk_l = 'Copy-Number genotyping from single cells ATAC/RNA'
    www = "https://github.com/caravagnalab/rcongas"
    em = "gcaravagna@gmail.com"

    cli::cli_alert_success(
      'Loading {.field {pk}}, {.emph \'{pk_l}\'}. Support : {.url { www}}' )

    options(rcongas_welcome_message = FALSE)
  }

  invisible()
}
