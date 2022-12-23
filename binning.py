import numpy as np

# Total number of binning steps
numbins = 13

# Read data from sample file
data = []
with open("thermodynamics.dat") as f:
    for line in f:
        step,Pot,Kin,Etot,Temp,Pre = map(float, line.strip().split(","))
        data.append({"x": Pre})

size0 = len(data)

# Initialize variables to compute averages 
x_avg = 0
averages = []

with open("binned.dat", "w") as f:
    for m in range(numbins):
        # Compute averages every 2**m steps
        bin_size = 2**m
        binned_data = [data[i:i+bin_size] for i in range(0, size0, bin_size)]
        for binned in binned_data:
                x_avg = sum(d["x"] for d in binned) / bin_size
                averages.append({"x": x_avg})

        # Compute averages and variance
        x_mean = np.mean([d["x"] for d in averages])
        x_std = np.std([d["x"] for d in averages])/np.sqrt(len(averages)-1)

        f.write(f"size = {2**m:5d} , x = {x_mean:7.5f} +- {x_std:7.5f}\n")
        averages = []
