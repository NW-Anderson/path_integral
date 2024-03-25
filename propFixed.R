setwd("/media/nathan/T7/path_integral/simulations/out/untarredFwdpySims")
prop_fixed <- c()
prop_lost <- c()
count <-0 
for(file in list.files()){
  count <- count + 1
  print(paste(count, file))
  start_freqs <- unlist(fread(paste(file,"/start_freqs.csv", sep=""), skip = 2)[,1])
  nsites <- length(start_freqs)
  nfixed <- sum(start_freqs == 1)
  nseg <- sum(start_freqs < 1 & start_freqs > 0)
  prop_fixed <- c(prop_fixed, nfixed / nseg)
  # prop_lost <- c(prop_lost, nlost / nsites)
}
