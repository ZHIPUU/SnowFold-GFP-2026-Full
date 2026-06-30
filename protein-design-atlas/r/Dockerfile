FROM rocker/tidyverse:4.4
RUN R -e "install.packages(c('DBI','RSQLite','patchwork','svglite','Cairo','plotly','htmlwidgets','igraph','tidygraph','ggraph','ComplexHeatmap','pheatmap','ggrepel'), repos='https://cloud.r-project.org')"
WORKDIR /app
COPY r /app/r
CMD ["Rscript", "/app/r/scripts/render_all.R"]
