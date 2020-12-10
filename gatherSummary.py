#Author: YP
#Created: 2018-08-24
#A scripts that gathers all of the results from the different scenarios
#into a single file, and then generates a summary of it.

#NOTE: This script hard-codes the ordering and contents of the csv file 
#summary. If you modify the code to change this, you will need to change
#it here too.

tsp = {'slices-tsp-rue-1000-3000':3600}
sat = {'slices-circuit-fuzz':600,'slices-BMC08':600}
mip = {'slices-cls':10000,'slices-regions200':10000}

scenarios = {'scenarios/LKH':tsp,'scenarios/EAX':tsp,'scenarios/CaDiCaL':sat,'scenarios/lingeling':sat,'scenarios/cryptominisat':sat,'scenarios/CPLEX':mip}

order = ['Parameter', 'Interesting', 'Unimodal', 'Convex', 'Trivially Convex', '% Instances Interesting', '% Instances Unimodal', '% Instances Convex', '% Instances Trivially Convex', 'Average number of modes per instance', 'Fitness Distance Correlation (Median)', 'Fitness Distance Correlation (Q2.5)', 'Fitness Distance Correlation (Q97.5)', 'Average Instance Fitness Distance Correlation (Median)', 'Average Instance Fitness Distance Correlation (Q2.5)', 'Average Instance Fitness Distance Correlation (Q97.5)']

with open('Parameter-Statistics-All-new.csv', 'w') as f_out:
  #data summary
  stats = []
  interestingStats = []
  #Write the headers
  f_out.write("#[1] Location")
  i = 2
  for term in order:
    f_out.write(', [' + str(i) + '] ' + term)
    i += 1
    stats.append([])
    interestingStats.append([])
  f_out.write('\n')

  for scenario in scenarios.keys():
    for phase in scenarios[scenario].keys():
      cutoff = scenarios[scenario][phase]
    
      with open(scenario + '/' + phase + '/responses/Parameter-Statistics.csv') as f_in:
        for line in f_in:
          if('#' in line[0]):
            continue
          f_out.write(scenario + '/' + phase + ', ' + line.strip() + '\n')
          items = line.split(',')

          interesting = items[1].strip() == '1'

          for i in range(0,len(items)):  
            try:
              items[i] = int(items[i])
            except:
              try:
                items[i] = float(items[i])
              except:
                pass
            stats[i].append(items[i])
            if(interesting):
              interestingStats[i].append(items[i])

with open('Parameter-Statistics-Summary-new.txt', 'w') as f_out:
  f_out.write('#Statistic: all (interesting)\n')
  for i in range(0,len(order)):
    if(type(stats[i][0]) is int):
      f_out.write('% ' + order[i] + ': ' + str(100.0*sum(stats[i])/len(stats[i])) +  ' (' + str(100.0*sum(interestingStats[i])/len(interestingStats[i])) + ')\n')
    elif(type(stats[i][0]) is float):
      f_out.write('Average ' + order[i] + ': ' + str(sum(stats[i])/len(stats[i])) + ' (' + str(sum(interestingStats[i])/len(interestingStats[i])) + ')\n')
    else:
      f_out.write(order[i] + ' count: ' + str(len(stats[i])) + ' (' + str(len(interestingStats[i])) + ')\n')
      
