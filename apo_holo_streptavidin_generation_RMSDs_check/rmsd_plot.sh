awk 'FNR>1{print $2}' unique_backbone.txt > plot_rmsd
xmgrace plot_rmsd
