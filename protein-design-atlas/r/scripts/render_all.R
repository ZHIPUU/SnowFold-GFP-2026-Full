#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DBI); library(RSQLite); library(ggplot2); library(dplyr); library(readr); library(patchwork)
})
db <- Sys.getenv('ATLAS_DB', '/app/data/design.db')
outdir <- Sys.getenv('ATLAS_PLOT_DIR', '/app/artifacts/plots')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
con <- dbConnect(SQLite(), db)
rounds <- dbGetQuery(con, "SELECT round_key, round_number, best_score FROM rounds WHERE best_score IS NOT NULL ORDER BY round_number")
seqs <- dbGetQuery(con, "SELECT source_round, best_score, best_ptm, best_plddt, best_chromo FROM sequences WHERE best_score IS NOT NULL")
dbDisconnect(con)

p1 <- ggplot(rounds, aes(round_number, best_score)) +
  geom_line(color='#d4af37', linewidth=1.2) + geom_point(color='#69c18d', size=3) +
  geom_text(aes(label=round_key), vjust=-1, size=3, color='#4b5563') +
  labs(title='Best score trajectory', x='Round', y='Best score') + theme_minimal(base_size=12)

p2 <- seqs %>% filter(best_score > 0.85) %>%
  ggplot(aes(best_ptm, best_chromo, color=best_score)) + geom_point(alpha=.55, size=1.2) +
  scale_color_viridis_c() + labs(title='pTM vs chromophore confidence', x='pTM', y='Chromo pLDDT') + theme_minimal(base_size=12)

p3 <- seqs %>% filter(best_score > 0.85) %>%
  ggplot(aes(source_round, best_score, fill=source_round)) + geom_boxplot(outlier.alpha=.15) +
  coord_flip() + guides(fill='none') + labs(title='Score distribution by round', x=NULL, y='Score') + theme_minimal(base_size=12)

combined <- (p1 / (p2 | p3)) + plot_annotation(title='Protein Design Atlas — R summaries')
ggsave(file.path(outdir, 'r_summary_dashboard.png'), combined, width=14, height=10, dpi=220)
ggsave(file.path(outdir, 'r_summary_dashboard.svg'), combined, width=14, height=10)
cat('wrote', file.path(outdir, 'r_summary_dashboard.png'), '\n')
