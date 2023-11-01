import fwdpy11
import numpy as np
from dataclasses import dataclass
from typing import List

import time
import demes

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
shift = 0.5

# Set up selected mutations
a = 0.01
sregions = [fwdpy11.ConstantS(0, L, 1, a), fwdpy11.ConstantS(0, L, 1, -a)]

# Load demographic model using demes
g = demes.load("subpopulations.yaml")
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
GSSmo = fwdpy11.GSSmo(
    [
        fwdpy11.Optimum(when=0, optimum=optimum, VS=VS),
        fwdpy11.Optimum(when=model.metadata["burnin_time"], optimum=shift, VS=VS),
    ]
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
        if pop.generation >= pop.N * burnin:
            sampler.assign(np.arange(0, pop.N))


recorder = Recorder(data=[])

## Initialize and evolve full population
pop = fwdpy11.DiploidPopulation(N0, L)
rng = fwdpy11.GSLrng(424242)

time1 = time.time()
fwdpy11.evolvets(rng, pop, params, 100, recorder=recorder, suppress_table_indexing=True)
time2 = time.time()
print(f"simulation took {(time2 - time1)/60:0.2f} minutes")
assert pop.generation == simlen

ts = pop.dump_tables_to_tskit()
line_time = g.demes[0].end_time
assert ts.num_samples == pop.N * 2 * (line_time + 1)


"""
## The tree sequence should record all generations from the last generation
## prior to the split into replicate lines, keeping every individual in all
## generations after the optimum shift.

print("Tree sequence has", ts.num_mutations, "mutations")
mut_times = []
for m in ts.mutations():
    mut_times.append(m.time)

# infinite sites, so this should be true
assert ts.num_mutations == ts.num_sites

# get times to compute frequencies at
times = set()
for s in ts.samples():
    times.add(ts.node(s).time)

times = sorted(list(times))[::-1]

allele_frequencies = {i: np.zeros((ts.num_sites, len(times))) for i in range(10)}

for j, t in enumerate(times):
    print(j)
    samples = [s for s in ts.samples() if ts.node(s).time == t]
    ts_slice = ts.simplify(samples, filter_sites=False)
    G = ts_slice.genotype_matrix()
    if t == max(times):
        # initial generation before shift in lines
        afs = G.sum(axis=1) / 2 / N0
        for i in range(10):
            allele_frequencies[i][:, j] = afs
    else:
        for i in range(10):
            G_sub = G[:, i * 2 * Nf : (i + 1) * 2 * Nf]
            afs = G_sub.sum(axis=1) / 2 / Nf
            allele_frequencies[i][:, j] = afs

#
"""
