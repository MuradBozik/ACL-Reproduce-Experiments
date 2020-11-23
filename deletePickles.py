import os

tsp = {'slices-tsp-rue-1000-3000':3600}
sat = {'slices-circuit-fuzz':600,'slices-BMC08':600}
mip = {'slices-cls':10000,'slices-regions200':10000}

scenarios = {'scenarios/EAX':tsp,'scenarios/LKH':tsp,'scenarios/CaDiCaL':sat,'scenarios/lingeling':sat,'scenarios/cryptominisat':sat,'scenarios/CPLEX':mip}


for scenario in scenarios.keys():
  for phase in scenarios[scenario].keys():
    cutoff = scenarios[scenario][phase]

    os.system('rm ' + scenario + '/' + phase + '/data.pkl -f')
    os.system('rm ' + scenario + '/' + phase + '/slice-data -rf')

