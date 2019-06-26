import cluster
import plots
import logging
import time


logging.getLogger().setLevel(logging.INFO)
name_folder = '/home/andrei/Documents/University/Computerpraktikum/Python/ExtractedFiles/Cluster/'
name = 'cod-rna.5000.csv'
output_folder = '/home/andrei/Documents/University/Computerpraktikum/Python/ExtractedFiles/Cluster/Artificial/Output/'
output_file = 'svmguide1.output'
time_file = open("/home/andrei/Documents/University/Computerpraktikum/Python/ExtractedFiles/Cluster/Artificial/Output/svmguide1.times.txt", "w")


for delta in [0.2, 0.1, 0.06, 0.05, 0.04]:
    for tau in [2*delta, 4*delta]:
        for eps_mul in [1, 1.5, 2]:
            ending = '._delta='+str(delta)+'_tau='+str(tau)+'_e='+str(eps_mul)+'_'
            stime = time.time()
            eps, cp = cluster.cluster(name=name_folder+name, delta=delta, epsilon=None, tau=tau,
                            return_list=True, epsilon_multiple=eps_mul, output_name=output_folder+output_file+ending)
            stime = time.time()-stime
            time_file.write("Time for "+output_file+ending+" with epsilon = " + str(eps) + " is: "+str(stime)+"\n")
            #plots.plot_components(cp, file=output_folder+output_file+ending+'.plot.png')

time_file.close()