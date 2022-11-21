analysis_name <- "k15_cov_enrichment"
source("/omics/groups/OE0219/internal/yuyu/10.Benchmarking/analysis/data_processing_redo/00.project_setting.R")
library(DescTools)

## read data:
################################################
# cummulative_cov_rank_tb <- read_("cummulative_cov_rank_tb", analysis_name)
# cummulative_cov_rank_condensed_tb <- read_("cummulative_cov_rank_condensed_tb", analysis_name)
# 
# beta_cov_tb <- read_("beta_cov_tb", "2_import_to_methrix") #%>% filter(chr=="chr1")
#     #filter(!(protocol == "EMseq" & sample == "6T")), 

beta_cov_tb <- read_(
  "cov_tb",
  "k04_meth_and_cov_distr"
)

## calc percent_rank:
########################

cov_rank_tb <- beta_cov_tb %>%
    dplyr::select(protocol, sample, chr, pos, cov) %>%
    group_by(sample, protocol) %>%
    dplyr::mutate(top_percent_rank=1-percent_rank(cov)) %>%
    ungroup
    
cov_rank_condensed_tb <- beta_cov_tb %>%
    dplyr::select(protocol, chr, pos, cov) %>%
    group_by(protocol) %>%
    dplyr::mutate(top_percent_rank=1-percent_rank(cov)) %>%
    ungroup

save_("cov_rank_condensed_tb", data=cov_rank_condensed_tb)

## cummulative coverage for given top_percent_rank:
###################################################
ranks_to_test <- seq(0,1,length.out=100)

total_counts_tb <- cov_rank_tb %>%
    group_by(protocol, sample) %>%
    summarize(total_count=sum(cov))

cummulative_cov_rank_tb <- lapply(ranks_to_test, function(cutoff){
    cov_rank_tb %>%
        filter(top_percent_rank <= cutoff) %>%
        group_by(protocol, sample) %>%
        summarize(count_above_rank=sum(cov)) %>%
        dplyr::mutate(top_percent_rank=cutoff)
}) %>%
    do.call(bind_rows, .) %>%
    left_join(total_counts_tb, by=c("protocol", "sample")) %>%
    dplyr::mutate(cummulative_fraction_of_counts=count_above_rank/total_count) %>%
    dplyr::select(protocol, sample, top_percent_rank, cummulative_fraction_of_counts) %>%
    set_order()

save_("cummulative_cov_rank_tb", data=cummulative_cov_rank_tb)


total_counts_condensed_tb <- cov_rank_condensed_tb %>%
    group_by(protocol) %>%
    summarize(total_count=sum(cov))

cummulative_cov_rank_condensed_tb <- lapply(ranks_to_test, function(cutoff){
    cov_rank_condensed_tb %>%
        filter(top_percent_rank <= cutoff) %>%
        group_by(protocol) %>%
        summarize(count_above_rank=sum(cov)) %>%
        dplyr::mutate(top_percent_rank=cutoff)
}) %>%
    do.call(bind_rows, .) %>%
    left_join(total_counts_condensed_tb, by=c("protocol")) %>%
    dplyr::mutate(cummulative_fraction_of_counts=count_above_rank/total_count) %>%
    dplyr::select(protocol, top_percent_rank, cummulative_fraction_of_counts) %>%
    set_protocol_order()


# set_protocol_order <- function(tb){
#   tb %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(protocol=ordered(protocol, levels=protocols))
# }
# 
# set_sample_order <- function(tb){
#   tb %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(sample=ordered(sample, levels=samples))
# }
# 
# set_order <- function(tb){
#   tb %>%
#     set_protocol_order() %>%
#     set_sample_order()
# }
# 
# a = cummulative_cov_rank_condensed_tb %>%
#   left_join(total_counts_condensed_tb, by=c("protocol")) %>%
#   dplyr::mutate(cummulative_fraction_of_counts=count_above_rank/total_count) %>%
#   select(protocol, top_percent_rank, cummulative_fraction_of_counts) %>%
#   #set_protocol_order()  %>%
#   
# 
#   set_sample_order()


save_("cummulative_cov_rank_condensed_tb", data=cummulative_cov_rank_condensed_tb)

# protocol_color_sheme <- c(
#   "#CD534CFF", #EB534C
#   "#33cc00",
#   "#0073C2FF", #0073C2
#   "#EFC000FF", #EFC000
#   "#ff9900"
# )
# names(protocol_color_sheme) <- protocols


## plot percent rank vs. cummulative fraction:
############################################################
custom_breaks <- seq(0.00, 1.00, 0.25)
custom_labels <- c("1.00", "0.75", "0.50", "0.25", "0.00")
top_percent_rank_vs_cum_frac_plot <- cummulative_cov_rank_tb %>%
    ggplot() +
        geom_line(
            aes(top_percent_rank, cummulative_fraction_of_counts, color=protocol),
            size=1,
            alpha=2
        ) +
        xlab("coverage rank") +
        ylab("cummulative fraction of counts") +
        scale_color_manual(values=protocol_color_sheme) +
        geom_abline(intercept=0, slope=1, size=1, linetype="dashed") +
        scale_x_continuous(breaks=custom_breaks, labels=custom_labels) +
        theme(legend.position="none", 
              panel.background = element_rect(fill = "white", colour = "black", linetype="solid"),
              strip.background =element_blank(),
              panel.border=element_rect(fill=NA),
              axis.text.x = element_text(angle = 45, hjust = 1)) +
        facet_grid(. ~ sample)

save_(
    "top_percent_rank_vs_cum_frac_plot",
    plot=top_percent_rank_vs_cum_frac_plot,
    width=7.5,
    height=3
)

top_percent_rank_vs_cum_frac_condensed_plot <- cummulative_cov_rank_condensed_tb %>%
    ggplot() +
        geom_line(
            aes(top_percent_rank, cummulative_fraction_of_counts, color=protocol),
            size=1,
            alpha=2
        ) +
        xlab("coverage rank") +
        ylab("cummulative fraction of counts") +
        scale_color_manual(values=protocol_color_sheme) +
        geom_abline(intercept=0, slope=1, size=1, linetype="dashed") +
        scale_x_continuous(breaks=custom_breaks, labels=custom_labels) +
        theme(panel.background = element_rect(fill = "white", colour = "black", linetype="solid"),
        panel.border=element_rect(fill=NA),
        strip.background = element_blank(),
        legend.position='none') 

top_percent_rank_vs_cum_frac_condensed_plot

save_(
    "top_percent_rank_vs_cum_frac_condensed_plot",
    plot=top_percent_rank_vs_cum_frac_condensed_plot,
    width=3.2,
    height=3.2
)

