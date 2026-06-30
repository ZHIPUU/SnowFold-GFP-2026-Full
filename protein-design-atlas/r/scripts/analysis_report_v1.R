#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})
root <- '/work'
data_dir <- file.path(root, 'artifacts', 'report_v1_data')
out_dir <- file.path(root, 'artifacts', 'report_v1')
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

rounds <- read.csv(file.path(data_dir, 'rounds.csv'), stringsAsFactors=FALSE)
seqs <- read.csv(file.path(data_dir, 'seq_scores.csv'), stringsAsFactors=FALSE)
top <- read.csv(file.path(data_dir, 'top_candidates.csv'), stringsAsFactors=FALSE)
ham <- read.csv(file.path(data_dir, 'hamming_top40.csv'), stringsAsFactors=FALSE)
parents <- read.csv(file.path(data_dir, 'parent_contrib.csv'), stringsAsFactors=FALSE)

# Fix types
for (nm in c('best_score','score','ptm','plddt','chromo','identity','hamming','max_score','mean_score','mean_ptm','mean_chromo')) {
  for (dfname in c('rounds','seqs','top','ham','parents')) {
    if (exists(dfname)) {
      df <- get(dfname); if (nm %in% names(df)) { df[[nm]] <- as.numeric(df[[nm]]); assign(dfname, df) }
    }
  }
}

classic <- theme_minimal(base_size=13) + theme(
  plot.title=element_text(face='bold', size=17, color='#1b2a24'),
  plot.subtitle=element_text(color='#5f6b63'),
  panel.grid.minor=element_blank(),
  axis.title=element_text(color='#283832'),
  axis.text=element_text(color='#35443e')
)

# 1. score trajectory
p1 <- ggplot(rounds, aes(round_number, best_score)) +
  geom_line(color='#254b45', linewidth=1.2) +
  geom_point(aes(size=best_score), color='#d4af37') +
  geom_text(aes(label=round_key), nudge_y=0.0025, size=3, color='#1c302a') +
  scale_y_continuous(limits=c(min(rounds$best_score, na.rm=TRUE)-0.01, max(rounds$best_score, na.rm=TRUE)+0.006)) +
  guides(size='none') + labs(title='Best score trajectory across rounds', subtitle='R22/R24/R26 form the current high-score frontier', x='Round', y='Best score') + classic

ggsave(file.path(out_dir,'01_score_trajectory.png'), p1, width=11, height=6.2, dpi=220)

# 2. scatter frontier
p2 <- ggplot(seqs %>% filter(!is.na(ptm), !is.na(chromo), score > 0.85), aes(ptm, chromo)) +
  geom_point(aes(color=score, size=score), alpha=.62) +
  geom_point(data=top[1:12,], aes(ptm, chromo), shape=21, fill=NA, color='#111111', size=3.5, stroke=.8) +
  scale_color_viridis_c(option='C', limits=c(0.85, max(seqs$score, na.rm=TRUE))) +
  scale_x_continuous(limits=c(0.75, 0.94)) + scale_y_continuous(limits=c(0.75, 0.98)) +
  labs(title='pTM × Chromophore pLDDT frontier', subtitle='Confidence axes normalized to 0–1; high-score cluster concentrates in upper-right', x='pTM', y='Chromophore pLDDT', color='Score') + classic

ggsave(file.path(out_dir,'02_ptm_chromo_frontier.png'), p2, width=10.5, height=7, dpi=220)

# 3. round distributions for high-score rounds
sel_rounds <- c('R19','R20','R22','R23','R24','R26')
p3 <- ggplot(seqs %>% filter(source_round %in% sel_rounds, score > 0.88), aes(source_round, score, fill=source_round)) +
  geom_violin(alpha=.55, color=NA) + geom_boxplot(width=.16, outlier.alpha=.08, color='#253a34') +
  geom_jitter(width=.1, alpha=.12, size=.8) + guides(fill='none') +
  labs(title='Score distribution by late-round strategy', subtitle='R24/R26 shift the upper envelope beyond R22', x=NULL, y='Score') + classic

ggsave(file.path(out_dir,'03_round_distribution.png'), p3, width=11, height=6.5, dpi=220)

# 4. top candidate component heatmap
heat <- top[1:30,]
heat$label <- paste0(heat$source_round, '_', seq_len(nrow(heat)), '\n', sprintf('%.4f', heat$score))
heat_long <- rbind(
  data.frame(label=heat$label, metric='pTM', value=heat$ptm),
  data.frame(label=heat$label, metric='pLDDT', value=heat$plddt),
  data.frame(label=heat$label, metric='Chromo', value=heat$chromo),
  data.frame(label=heat$label, metric='Score', value=heat$score)
)
heat_long$label <- factor(heat_long$label, levels=rev(heat$label))
p4 <- ggplot(heat_long, aes(metric, label, fill=value)) + geom_tile(color='white', linewidth=.25) +
  scale_fill_viridis_c(option='A') + labs(title='Top 30 candidate metric heatmap', subtitle='R26/R24 candidates dominate the current score ceiling', x=NULL, y=NULL, fill='value') + classic + theme(axis.text.y=element_text(size=7))

ggsave(file.path(out_dir,'04_top30_heatmap.png'), p4, width=8.5, height=10, dpi=220)

# 5. hamming identity matrix
ham$seq_i <- factor(ham$seq_i, levels=unique(ham$seq_i))
ham$seq_j <- factor(ham$seq_j, levels=rev(unique(ham$seq_j)))
p5 <- ggplot(ham, aes(seq_i, seq_j, fill=identity)) + geom_tile() +
  scale_fill_viridis_c(option='B', limits=c(0,1)) + labs(title='Top 40 sequence identity matrix', subtitle='Block structure highlights local exploitation around successful parents', x=NULL, y=NULL, fill='identity') + classic + theme(axis.text.x=element_text(angle=90, hjust=1, size=6), axis.text.y=element_text(size=6))

ggsave(file.path(out_dir,'05_top40_identity_matrix.png'), p5, width=10, height=9, dpi=220)

# 6. parent contribution
parents2 <- parents %>% filter(!is.na(max_score)) %>% arrange(desc(max_score)) %>% head(18)
parents2$parent_label <- paste0(parents2$source_round, ' / ', substr(parents2$parent,1,18))
parents2$parent_label <- factor(parents2$parent_label, levels=rev(parents2$parent_label))
p6 <- ggplot(parents2, aes(max_score, parent_label)) +
  geom_segment(aes(x=mean_score, xend=max_score, y=parent_label, yend=parent_label), color='#8aa39b', linewidth=2.2, alpha=.8) +
  geom_point(aes(size=n, color=mean_chromo), alpha=.95) +
  scale_color_viridis_c(option='D') + labs(title='Parent contribution: mean → max score', subtitle='Best parents both improve score and preserve chromophore confidence', x='Score', y=NULL, size='n', color='mean chromo') + classic

ggsave(file.path(out_dir,'06_parent_contribution.png'), p6, width=11, height=7.5, dpi=220)

# composite
combo <- (p1 | p2) / (p3 | p6)
ggsave(file.path(out_dir,'00_composite_dashboard.png'), combo, width=17, height=12, dpi=190)

# HTML report
best <- top[1,]
imgs <- c('00_composite_dashboard.png','01_score_trajectory.png','02_ptm_chromo_frontier.png','03_round_distribution.png','04_top30_heatmap.png','05_top40_identity_matrix.png','06_parent_contribution.png')
img_tags <- paste0('<section class="figure"><h2>', seq_along(imgs), '. ', imgs, '</h2><img src="report_v1/', imgs, '"/></section>', collapse='\n')
html <- paste0('<!doctype html><html><head><meta charset="utf-8"><title>Protein Design Atlas — Analysis Report v1</title><style>',
'body{margin:0;background:#0e1110;color:#efe8d8;font-family:Georgia,serif} header{padding:56px 70px;background:linear-gradient(135deg,#24332d,#0e1110);border-bottom:1px solid #3b4a42} h1{font-size:54px;margin:0;color:#d4af37} .sub{color:#b9b09c;font-family:monospace}.grid{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;padding:24px 70px}.card{background:#151b18;border:1px solid #34443e;border-radius:18px;padding:18px}.card b{display:block;font-size:32px;color:#d4af37}.content{padding:10px 70px 80px}.figure{background:#f7f4eb;color:#17211d;border-radius:22px;padding:18px;margin:26px 0}.figure img{width:100%;border-radius:14px;box-shadow:0 10px 30px #0003} code{color:#9fdbb2}',
'</style></head><body><header><div class="sub">R4–R26 / SQLite + R graphics / report_v1</div><h1>Protein Design Atlas<br/>Analysis Report v1</h1><p class="sub">Current best: ', best$source_round, ' score=', sprintf('%.4f', best$score), ' pTM=', sprintf('%.4f', best$ptm), ' chromo=', sprintf('%.4f', best$chromo), '</p></header>',
'<div class="grid"><div class="card">Sequences<b>', nrow(seqs), '</b></div><div class="card">Rounds<b>', nrow(rounds), '</b></div><div class="card">Top score<b>', sprintf('%.4f', best$score), '</b></div><div class="card">Top chromo<b>', sprintf('%.4f', max(seqs$chromo,na.rm=TRUE)), '</b></div></div>',
'<div class="content"><p>本报告由 R/ggplot2 自动生成，图表已基于统一后的 0–1 pLDDT/chromophore 尺度。重点显示 R22→R24→R26 的高分前沿、各轮分布、Top 候选热图、序列相似性矩阵与父代贡献。</p>', img_tags, '</div></body></html>')
writeLines(html, file.path(root,'artifacts','analysis_report_v1.html'))
cat('report:', file.path(root,'artifacts','analysis_report_v1.html'), '\n')
