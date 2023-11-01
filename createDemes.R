lines <- c("time_units: generations",
           "defaults:",
           "  epoch: {start_size: 200, end_time: 0}",
           "demes:",
           "  - name: ancestral",
           "epochs:",
           "  - {start_size: 2000, end_time: 100}",
           paste("- {name: line",
                 1:1000,
                 ", ancestors: [ancestral]}", 
                 sep = ""))

fileConn<-file("demo.yaml")
writeLines(lines, fileConn)
close(fileConn)

