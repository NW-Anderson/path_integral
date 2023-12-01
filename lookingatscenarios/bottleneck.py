import fwdpy11
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List
import random
import csv
import matplotlib.pyplot as plt



import time
import demes

import os

for seed in random.sample(range(1,int(1e5)), 3):
    ## Set up parameters
    
    # 1 Morgan chrosome
    r = 1e-8
    L = 1e8
    assert r * L == 1
    
    # VS always set to 1
    VS = 1.0
    
    # Change U in order to adjust VG
    U = 2.5e-3
    
    expectedVG = 4 * U * VS
    optimum = 0.0
    shift = 1.0
    
    # Set up selected mutations
    a = 0.01
    sregions = [fwdpy11.ConstantS(0, L, 1, a), fwdpy11.ConstantS(0, L, 1, -a)]
    
    # Load demographic model using demes
    os.chdir("/home/nathan/Documents/GitHub/path_integral")
    g = demes.load("subpopulations.yaml")
    
    # import demesdraw
    
    # demesdraw.tubes(g, colours=("blue"))
    burnin = 10
    model = fwdpy11.discrete_demography.from_demes(g, burnin=burnin)
    simlen = model.metadata["total_simulation_length"]
    
    # Setting up the population
    initial_sizes = [
        model.metadata["initial_sizes"][i]
        for i in sorted(model.metadata["initial_sizes"].keys())
    ]
    N0 = initial_sizes[0]
    assert len(initial_sizes) == 1
    Nf = g.demes[-1].epochs[0].start_size
    
    pop = fwdpy11.DiploidPopulation(initial_sizes, L)
    
    # after the burnin, the population is split into
    # the replicate lines and the optimum shift occurs
    GSSmo = fwdpy11.GSS(
          fwdpy11.Optimum(when=0, optimum=optimum, VS=VS)
    )
    
    pdict = {
        "nregions": [],
        "sregions": sregions,
        "recregions": [fwdpy11.BinomialInterval(0, L, 1)],
        "rates": (0.0, U, None),
        "gvalue": fwdpy11.Additive(scaling=2, gvalue_to_fitness=GSSmo),
        "simlen": simlen,
        "demography": model,
        "prune_selected": False,
    }
    params = fwdpy11.ModelParams(**pdict)
    
    ## set up recorders
    @dataclass
    class SimData:
        generation: int
        demes_ids: List[int]
        mean_phenotype: List[float]
        mean_fitness: List[float]
        var_phenotype: List[float]
    
    
    @dataclass
    class Recorder:
        data: list
    
        def __call__(self, pop, sampler):
            md = np.array(pop.diploid_metadata)
            # general properties of the population
            # store lists of mean phenotypes and fitness, and var(pheno)
            deme_ids = sorted(list(set(md["deme"])))
            mean_pheno = [md[md["deme"] == i]["g"].mean() for i in deme_ids]
            mean_fitness = [md[md["deme"] == i]["w"].mean() for i in deme_ids]
            var_pheno = [md[md["deme"] == i]["g"].var() for i in deme_ids]
            self.data.append(
                SimData(pop.generation, deme_ids, mean_pheno, mean_fitness, var_pheno)
            )
            # record last generation of full population and
            # every generation of the replicate lines
            if pop.generation == pop.N * burnin or pop.generation == simlen:
                sampler.assign(np.arange(0, pop.N))
    
    
    
    recorder = Recorder(data=[])
    print(seed)
    ## Initialize and evolve full population
    pop = fwdpy11.DiploidPopulation(N0, L)
    rng = fwdpy11.GSLrng(seed)

    time1 = time.time()
    fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder, suppress_table_indexing=True)
    time2 = time.time()
    print(f"simulation took {(time2 - time1)/60:0.2f} minutes")
    assert pop.generation == simlen

    ts = pop.dump_tables_to_tskit()
    line_time = g.demes[0].end_time
    
    gen = []
    mean_pheno = []
    mean_fit = []
    var_pheno = []
    for i in range(simlen): 
        gen.append(float(recorder.data[i].generation))
        mean_pheno.append(recorder.data[i].mean_phenotype[0])
        mean_fit.append(recorder.data[i].mean_fitness[0])
        var_pheno.append(recorder.data[i].var_phenotype[0])
    
    initial_VG = np.mean(var_pheno[N0 * (burnin - 1):N0 * burnin])
    
    f, ax = plt.subplots()
    ax.plot(gen, mean_pheno, label="Mean Phenotype")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Value")
    plt.legend()
    plt.show()
    
    print(initial_VG)
    print(f"The tree sequence now has {ts.num_mutations} mutations, at "
          f"{ts.num_sites} distinct sites.")
    
    with open(str(seed) + "popStats.csv", "w") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["scenario","gen","mean_pheno",  "var_pheno","mean_fit"])
        for i in range(len(gen)):
            writer.writerow(["bottleneck", gen[i], mean_pheno[i], var_pheno[i], mean_fit[i]])


    
    
