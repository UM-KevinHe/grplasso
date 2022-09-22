setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/Runtime")

plt.RMSE <- annotate_figure(ggarrange(Runtime.plot.pplasso,
                                      Runtime.plot.grplasso,
                                      nrow = 1, ncol = 2, 
                                      labels = c("A", "B")),
                            bottom = text_grob("Number of providers", face = "bold", 
                                               family = "serif", size = 14),
                            left = text_grob("Runtime (mins)", face = "bold", 
                                             family = "serif", size = 14, rot = 90))

save(plt.RMSE,
     file = paste0("Runtime_combine_", Sys.Date(), ".RData"))