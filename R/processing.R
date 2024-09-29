#' theme_border
#'
#' ggplot2 aesthetic option with borders
#' @export
theme_border <- function(){
    theme(
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
        }

#' theme_line
#'
#' ggplot2 aesthetic option with borders
#' @export
theme_line <- function(){
    list(
        scale_y_continuous(expand = c(0, 0)),
        theme(
            axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1),
            panel.grid = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", size = 1)))
        }

#' facet_aes
#'
#' ggplot2 aesthetic option with facet_wraps
#' @export
facet_aes <- function(){
    theme(
	strip.background = element_blank(),
        strip.text = element_text(face="bold", size=12))}


#' collapse_reads
#' 
#' Collapse shallow/deep sequencing reads
#' 
#' @param x dataframe of JH4 pipeline output
#' @export
collapse_reads <- function(x){
    x <- x %>%
        mutate(
            indel = case_when(
                deletions + insertions > 0 ~ TRUE,   # if there are any deltions and insertions
                .default = FALSE)) %>%
            group_by(sample, snps, deletions, insertions, total, length, sequence, conversions, indel) %>%
            summarize(count = sum(count)) %>%   # collapse reads from same library (shallow/deep sequencing)
            ungroup() %>%
            arrange(desc(count)) %>%
            mutate(rank = dense_rank(dplyr::desc(count))) %>%   # calculate new rank
            as.data.frame(.)
    return(x)}


#' plot_seq_length
#' 
#' plot the sequencing length from each sample
#' 
#' @param x dataframe of JH4 pipeline output
#' @param sample column name indicating sample
#' @param group column name indicating group
#' @export
plot_seq_length <- function(x, sample = "sample", group = NULL, ...){
    plot <- x %>%
        ggplot(aes_string(x = paste0(sample), y="length", color = group)) +
        geom_jitter(size = 2, width = 0.2) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        geom_hline(yintercept = 433) +
        xlab("") +
        ylab("Sequence length (bp-1)") +
        theme(
            axis.text.x = element_text(angle = 20, vjust = 0.8, hjust=0.8),
            panel.grid = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.8))

    if(length(group) > 0){
        facet_var <- paste0('~ ', group)
        plot <- plot + facet_wrap(as.formula(facet_var), nrow = 1, scales = "free")}
    
    return(plot)}

#' remove_indels
#' 
#' remove sequences with insertions and deletions
#' 
#' @param x dataframe of JH4 pipeline output
#' @param sample column name indicating sample
#' @param filter whether to remove indel sequences in the final output
#' @export
remove_indels <- function(x, sample = "sample", filter = T){
    x %>%
        mutate(indel = ifelse(insertions + deletions > 0, "WITH_INDEL", "NO_INDEL")) %>%
            group_by_at(vars(one_of(c(sample, "indel")))) %>%
            summarize(count = n()) %>%
            group_by_at(vars(one_of(sample))) %>%
            mutate(pct = count*100/sum(count)) %>%
            dplyr::select(!count) %>%
            pivot_wider(names_from = "indel", values_from = "pct") %>%
            print(.)
    if(filter == T){
        x <- x %>%
            mutate(indel = ifelse(insertions + deletions > 0, "WITH_INDEL", "NO_INDEL")) %>%
            filter(indel == "NO_INDEL")}
    return(x)}

#' plot_multi_qc
#' 
#' plot quality control parameters
#' 
#' @param x dataframe of JH4 pipeline output
#' @param sample column name indicating sample
#' @param filter whether to remove indel sequences in the final output
#' @export
plot_multi_qc <- function(x, sample = "sample", group = NULL){
    x <- x %>%
        group_by_at(vars(one_of(sample, group))) %>%
        summarize(nseq = n(), nreads = sum(count), nsnps = sum(snps), msnps = mean(snps))

    size <- 4
    if(length(group) > 0){
        color <- group}
    else{
        color <- var}

    p1 <- x %>%
        ggplot(aes(x=nreads, y=nseq)) +
            geom_point(aes_string(color=color), size = size) +
            ylab("nSeq") +
            xlab("nReads") +
            geom_vline(xintercept = 400, color = "red", linetype = "dashed") +
            theme_border() +
            theme(aspect.ratio = 1) +
            guides(color = guide_none(), shape = guide_none())

    p2 <- x %>%
        ggplot(aes(x=nreads, y=nsnps)) +
            geom_point(aes_string(color=color), size = size) +
            ylab("nSNPs") +
            xlab("nReads") +
            geom_vline(xintercept = 400, color = "red", linetype = "dashed") +
            theme_border() +
            theme(aspect.ratio = 1) +
            guides(color = guide_none(), shape = guide_none())

    p3 <- x %>%
        ggplot(aes(x=nseq, y=nsnps)) +
            geom_point(aes_string(color=color), size = size) +
            ylab("nSNPs") +
            xlab("nSeq") +
            theme_border() +
            theme(aspect.ratio = 1) +
            guides(color = guide_none(), shape = guide_none())


    p4 <- x %>%
        ggplot(aes(x=nseq, y=msnps)) +
            geom_point(aes_string(color=color), size = size) +
            ylab("AvgSNPs") +
            xlab("nSeq") +
            theme_border() +
            theme(aspect.ratio = 1) +
            guides(color = guide_none(), shape = guide_none())
    p5 <- x %>%
        ggplot(aes(x=nreads, y=msnps)) +
            geom_point(aes_string(color=color), size = size) +
            geom_vline(xintercept = 400, color = "red", linetype = "dashed") +
            ylab("AvgSNPs") +
            xlab("nReads") +
            theme_border() +
            theme(aspect.ratio = 1)

    plot <- plot_grid(p1, p2, p3, p4, p5, ncol = 3, align = "hv")
    return(plot)
}

#' plot_seq_copy
#' 
#' plot number of reads for each unique sequence
#' 
#' @param x dataframe of JH4 pipeline output
#' @param var column name to group by
#' @export
plot_seq_copy <- function(x, var){
    plot <- x %>%
        ggplot(aes(x=log10(count), y=total)) +
        geom_point(aes_string(colour=var)) +
        geom_vline(xintercept = log10(1.5), linetype = "dashed") +
        xlab("log10 nReads\nin each sequence") +
        ylab("nSNPs\nin each sequence") +
        theme_border()
    return(plot)}

#' remove_lowqc
#' 
#' remove samples that do not pass quality control
#' 
#' @param x dataframe of JH4 pipeline output
#' @param sample column name indicating sample
#' @param min.reads minimum number of reads; defaults to 400
#' @param min.copy minimum reads for each sequence; defaults to 1
#' @export
remove_lowqc <- function(x, sample = "sample", min.reads = 400, min.copy = 1){
    # pass min.reads
    x <- x %>%
        group_by_at(vars(one_of(c(sample)))) %>%
        mutate(pass.min.reads = ifelse(sum(count) >= min.reads, "Pass", "Fail")) %>%
        group_by_at(vars(one_of(c(sample, "sequence")))) %>%
        mutate(pass.min.copy = ifelse(count >= min.copy, "Pass", "Fail")) %>%
        filter(pass.min.reads == "Pass", pass.min.copy == "Pass") %>%
        ungroup()

    if(is.factor(x[[sample]])){
        levels <- levels(x[[sample]])[which(levels(x[[sample]]) %in% unique(x[[sample]]))]
        x[[sample]] <- factor(x[[sample]], levels)}

    return(x)}


#' count_snp
#' 
#' count number of SNP
#' 
#' @param x dataframe of JH4 pipeline output
#' @param absolute defaults to TRUE to return SNP counts from all reads, else return SNP counts from each sequence
#' @export
count_snp <- function(x, absolute = T){

    # make snp table
    table <- matrix(nrow=4, ncol=4, 0)
    colnames(table) <- c('A', 'T', 'C', 'G')
    rownames(table) <- c('A', 'T', 'C', 'G')

    # fill snp table
    conversions <- strsplit(x$conversions,';')

    for (j in 1:length(conversions)){
        nreads <- as.numeric(x$count[j])
        for (k in 1:length(conversions[[j]])){
            snp_count <- unlist(strsplit(conversions[[j]][k],'='))
            from_to <- unlist(strsplit(snp_count[1],':'))

            # if use absolute or not (read- or sequence- level)
            if(absolute == F){
                table[from_to[1], from_to[2]] <- table[from_to[1], from_to[2]] + as.numeric(snp_count[2])}
            else{
               	table[from_to[1], from_to[2]] <- table[from_to[1], from_to[2]] + as.numeric(snp_count[2])*nreads}}}

    # return snp table
    return(table)}

#' calculate_sample_snp
#' 
#' calculate SNP count matrix per sample
#' 
#' @param x dataframe of JH4 pipeline output
#' @param sample column name indicating sample
#' @param absolute defaults to TRUE to return SNP counts from all reads, else return SNP counts from unique sequences
#' @export
calculate_sample_snp <- function(x, sample = "sample", absolute = T){
    snptable <- list()
    samplelist <- split(x, x[[sample]])

    # calculate snp per sample
    for(i in seq_along(samplelist)){
        snptable[[i]] <- count_snp(samplelist[[i]], absolute = absolute)}

    names(snptable) <- names(samplelist)
    return(snptable)}

#' calculate_group_snp
#' 
#' calculate SNP count matrix per group
#' 
#' @param snptable a list of count matrix per sample
#' @param group_vec a vector of groups for each sample; vector must be same length as the list of snptable
#' @param average whether to average the counts per sample; defaults to TRUE
#' @export
calculate_group_snp <- function(snptable, group_vec = NULL, average = T){
    if(length(group_vec) == 0){
        group_vec <- rep("group", length(snptable))}
    stopifnot(length(snptable) == length(group_vec))

    avgsnptable <- list()
    glevels <- unique(group_vec)
    for(g in seq_along(glevels)){
        grouptable <- snptable[which(group_vec == glevels[g])]
        avgsnptable[[g]] <- Reduce("+", grouptable)
        if(average){
            avgsnptable[[g]] <- avgsnptable[[g]]/length(grouptable)}}
    names(avgsnptable) <- glevels

    return(avgsnptable)
}


#' plot_snp_heatmap
#' 
#' plot heatmap for SNP count matrices
#' 
#' @param snp_matrix_list a list of SNP matrix
#' @export
plot_snp_heatmap <- function(snp_matrix_list){
    heatmap_list <- list()
    percent_list <- list()

    legend <- rep(F, length(snp_matrix_list))
    legend[1] <- T

    for(x in seq_along(snp_matrix_list)){
        snptable <- snp_matrix_list[[x]]
        percent <- snptable*100/sum(snptable)
        rownames(percent) <- rownames(snptable)
        percent_list[[x]] <- percent}

    max <- Reduce("max", percent_list)

    # make heatmap
    for(x in seq_along(snp_matrix_list)){
        snptable <- snp_matrix_list[[x]]
        percent <- percent_list[[x]]
        col_fun = circlize::colorRamp2(seq(from = 0, to = round(max, -1), length.out = 10), rev(RColorBrewer::brewer.pal(10, "RdBu")))
        percent[percent == 0] <- NA
        heatmap_list[[x]] <- Heatmap(
            percent,
            border = TRUE,
            column_title = names(snp_matrix_list)[x],
            na_col = "black",
            col = col_fun,
            rect_gp = gpar(col = "white", lwd = 1),
            heatmap_legend_param = list(
                title = "SNP Count (%)",
                legend_direction = "vertical",
                width = unit(10, "mm")),
            show_heatmap_legend = legend[x],
            show_column_names = T,
            column_names_side = "bottom",
            #column_title = "To",
            column_title_side = "top",
            column_order = colnames(snptable),
            column_dend_reorder = F,
            show_row_names = T,
            row_title = "From",
            row_order = colnames(snptable),
            row_dend_reorder = F,
            row_names_side = "left",
            column_names_rot = 0,
            column_title_gp = gpar(fontsize = 15, fontface = "bold"),
            column_names_gp = gpar(fontsize = 10),
            row_title_gp = gpar(fontsize = 12),
            row_names_gp = gpar(fontsize = 10),
            cell_fun = function(j, i, x, y, w, h, fill){
                grid.text(round(snptable[i, j]), x, y)})
        }

    names(heatmap_list) <- names(snp_matrix_list)
    formula <- paste0("heatmap_list[[", 1:length(heatmap_list), "]]", collapse = " + ")
    expression_string <- paste0('draw(', formula, ', heatmap_legend_side = "right")')
    return(eval(parse(text = expression_string)))}

#' count_transitions
#' 
#' count the number of SNPs for transition frequency
#' 
#' @param conversions 
#' @param count 
#' @param absolute defaults to TRUE to return SNP counts from all reads, else return SNP counts from unique sequences
#' @export
count_transitions <- function(conversions, count, absolute = T){
    conversions_list <- strsplit(conversions,';')
    transition_list <- list()

    # specify each SNP
    for(position in 1:length(conversions_list[[1]])){

        # count this SNP in every sequence
        snp <- c()
        for(i in 1:length(conversions_list)){
            snp <- c(snp, conversions_list[[i]][position])}

        # if absolute = T, multiply nSNPs by nReads in this sequence
        if(absolute == F){
            snp_count <- as.numeric(unlist(str_extract_all(snp, "\\d+\\.?\\d*")))}
        else{
            snp_count <- as.numeric(unlist(str_extract_all(snp, "\\d+\\.?\\d*")))*count}

        # calculate sum for this SNP
        transition_list[[position]] <- sum(snp_count)}

    names(transition_list) <- gsub("\\=.*", "", conversions_list[[1]])
    return(transition_list)}

#' calculate_group_transition_freq
#' 
#' calculate SNP transition frequency per group
#' 
#' @param x dataframe of JH4 pipeline output
#' @param var column name to group by
#' @param absolute defaults to TRUE to return SNP counts from all reads, else return SNP counts from unique sequences
#' @export
calculate_group_transition_freq <- function(x, var, absolute = T) {
  df <- x %>%
    group_by(across(all_of(c(var)))) %>%
    summarize(
      nseq = n(),
      nreads = sum(count),
      transition_counts = list(count_transitions(conversions, count, absolute = absolute)),
      .groups = 'drop') %>%
    unnest_longer(transition_counts, indices_to = "SNP")

  if(absolute == F){
    df <- df %>%
      mutate(transition_freq = transition_counts/nseq)}
  else{
    df <- df %>%
      mutate(transition_freq = transition_counts/nreads)}

  return(df)
}

#' calculate_mu_freq
#' 
#' calculate SHM frequency
#' 
#' @param x dataframe of JH4 pipeline output
#' @param var column name to group by
#' @export
calculate_mu_freq <- function(x, var){
    df <- x %>%
        mutate(mu = total*count) %>%
        group_by(across(all_of(c(var)))) %>%
        summarize(nreads = sum(count), mu_freq = (sum(mu)/nreads)/433) %>%
        mutate(
            mu_freq_sf = mu_freq,
            mu_freq_sf = format(mu_freq_sf, scientific = TRUE, digits = 3),
            mu_freq_sf = paste(mu_freq_sf, " bp", "-1", sep = ""),
            mu_freq_sf = gsub("e", "x 10", mu_freq_sf),
            mu_freq_sf = gsub("-01", "\u207B\u00B1", mu_freq_sf, fixed = TRUE),
            mu_freq_sf = gsub("-02", "\u207B\u00B2", mu_freq_sf, fixed = TRUE),
            mu_freq_sf = gsub("-03", "\u207B\u00B3", mu_freq_sf, fixed = TRUE),
            mu_freq_sf = gsub("-1", "\u207B\u00B9", mu_freq_sf, fixed = TRUE)) %>%
        dplyr::select(c(var, "mu_freq", "mu_freq_sf"))

    return(df)}

#' plot_mu_prop
#' 
#' plot proportion of reads with SHM
#' 
#' @param x dataframe of JH4 pipeline output
#' @param var column name to group by
#' @export
plot_mu_prop <- function(x, var, ...){

    # split sequence into mu groups
    x$mu_groups <- ifelse(x$total >= 7, "≥7", x$total)
    x$mu_groups <- factor(x$mu_groups, c(0:6, "≥7"))

    # set colors for donut plot
    colors <- viridis(8)
    names(colors) <- levels(x$mu_groups)
    donuts <- list()

    # calculate mu_freq
    x$var_tmp <- x[[var]]
    mu_table <- calculate_mu_freq(x = x, var = "var_tmp")

    # calculate donut percentage
    shm <- x %>%
        group_by(var_tmp, mu_groups) %>%
        summarize(count = sum(count)) %>%
        group_by(var_tmp) %>%
        mutate(nreads = sum(count), pct = count/sum(count)) %>%
        mutate(ymax = cumsum(pct), ymin = c(0, head(ymax, n=-1)))

    # make donut ggplot
    if(is.factor(x$var_tmp)){
        varlevels <- levels(x$var_tmp)}
    else{
        varlevels <- unique(x$var_tmp)}

    for(i in seq_along(varlevels)){
        freq = mu_table %>% filter(var_tmp == varlevels[i]) %>% .$mu_freq_sf
        donuts[[i]] <- shm %>%
            filter(var_tmp == varlevels[i]) %>%
            ggplot(aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=mu_groups)) +
                geom_rect() +
                geom_text(x=2, aes(y=0, label=nreads), size=6) +
                scale_fill_manual(values = colors) +
                coord_polar(theta="y") +
                xlim(c(2, 4)) +
                theme_void() +
                theme(legend.position = "none") +
                ggtitle(paste0(varlevels[i]), subtitle = freq)}

    # plot legend
    legend <- ggplot(shm, aes(x = ymin, y = ymax, fill = mu_groups)) +
        geom_col() +
        scale_fill_manual(values = colors) +
        guides(fill = guide_legend(title = "No. of Mutations"))
    donuts[[length(donuts) + 1]] <- get_legend(legend)

    # merge
    plot <- plot_grid(plotlist = donuts, ...)
    return(plot)}