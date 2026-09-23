import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import subprocess, gzip
import numpy as np

base_pop_size = 1000000
q_val = 250

pop_size = base_pop_size / q_val
r_rate = 1e-8 * q_val
u_rate = 3e-9 * q_val
length = 10e+3

pn = 0.5
pb = 0.5 * 0.0002
pd = 0.5 * (1 - 0.0002)

bs = 250.0 / 4.0 / pop_size
ds = 100.0 / 4.0 / pop_size

t = 14 * pop_size

subprocess.run("../slim" + " -m" + " -l 0" +
               " -d" + " N=" + str(pop_size) +
               " -d" + " r=" + str(r_rate) +
               " -d" + " u=" + str(u_rate) +
               " -d" + " L=" + str(length) +
               " -d" + " pn=" + str(pn) +
               " -d" + " pb=" + str(pb) +
               " -d" + " pd=" + str(pd) +
               " -d" + " sb=" + str(bs) +
               " -d" + " sd=" + str(ds) +
               " -d" + " simend=" + str(t) +
               " simulator.slim", shell=True
)

output = gzip.open("fixed_mutation_list.txt.gz", mode="wt")

with open("substitutions.txt", mode="r") as f:
  for line in f:
    line = line.strip("\t")

    if line.startswith("#") or line.startswith("Mutations"):
      continue
    else:
      cols = line.split()
      print(cols[8], cols[7], cols[3], np.log(1.0 + float(cols[4])), file=output)

output.close()

exit()