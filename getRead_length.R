
real_data_read_length <- read.table('/local_scratch/real_data_read_length.txt',header = FALSE)

real_data_read_mapped_length <- read.table('/local_scratch/real_data_read_mapped_length.txt',header = FALSE)

hist(real_data_read_mapped_length$V1)

Simulation_data_read_length <- read.table('/local_scratch/Simulation_data_read_length.txt',header = FALSE)

Simulation_data_read_mapped_length <- read.table('/local_scratch/Simulation_data_read_mapped_length.txt',header = FALSE)

hist(Simulation_data_read_length$V1)

hist(Simulation_data_read_mapped_length$V1)