# Utiiliser une image docker existant déjà avec la version volue
FROM rocker/r-base:4.2.1

# On demande à installer les packages BiocManager (ddépendence de DESeq2)
#et DESeq2
Run R -e 'install.packages("BiocManager")'
Run R -e 'BiocManager::install("DESeq2")'
