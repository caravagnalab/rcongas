# 
# ggsave(
#   plot = report_analysis_chr(x),
#   filename = 'a.png',
#   width = 10,
#   height = 10)
  
#' Title
#'
#' @param x
#' @param ...
#' @param cex 
#'
#' @return
#' @export
#'
#' @importFrom ggpubr ggarrange
#' @importFrom CNAqc plot_segments
#' @importFrom patchwork plot_layout
#'
#' @examples
report_analysis_chr <- function(x, 
                                normalise = 'ls',
                                top = 5,
                                n_sorrounding = 5,
                                cex = 1, ...)
{
  highlights <-
    x %>% Rcongas::get_clusters_ploidy() %>% filter(highlight == TRUE)
  
  chr <- unique(highlights$chr)
  
  curate = function(p) {
    p +
      theme(
        plot.title = element_text(
          size = 11 * cex,
          color = "indianred3",
          face = "bold",
          family = "Roboto Black"
        ),
        plot.subtitle = element_text(color = "#666666", size = 7 * cex),
        plot.caption = element_text(color = "#AAAAAA", size = 5 * cex)
      )
  }
  
  # All plots
  list_plots = lapply(
    chr,
    plot_chromosome_inspection,
    normalise = normalise,
    top = top,
    n_sorrounding = n_sorrounding,
    x = x
  )
  
  
  # arrange
  # nplots = list_plots %>% length
  # 
  # nr = nplots %>% sqrt %>% ceiling()
  # nc = nplots %>% sqrt %>% ceiling()
  # 
  # if(nplots < 4)
  # {
  #   nr = 1
  #   nc = nplots
  # }
  
  nc = 1
  nr = list_plots %>% length
    
  
  ggpubr::ggarrange(
    plotlist = list_plots,
    nrow = nr,
    ncol = nc
  )
  
  
}