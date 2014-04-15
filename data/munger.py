import csv
import os  
import numpy as np

#fname = "raw/unc.edu.2d67443d-14ec-48f7-93e0-e409b1729587.1201146.rsem.genes.results"

def get_path(name, base=os.path.dirname(__file__)):
    return os.path.join(base, name)

kirc = get_path('raw')

counts = {}
genes = {}

for f in os.listdir(kirc):
     ff = get_path(f, kirc)
     if os.path.isfile(ff):
        subject = ff[-26:][:7]  ## this appears to be a unique identifier, but need to verify
        print subject
        counts[subject] = np.loadtxt(ff, skiprows=1, dtype=np.str, usecols=(1,))
        genes[subject] = np.loadtxt(ff, skiprows=1, dtype=np.str, usecols=(0,))




with open('data.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    subjects = counts.keys()
    writer.writerow(tuple(['gene']+subjects))
    data  = np.array(counts.values())
    for gene, count in zip(genes.values()[0], data.T):
        writer.writerow(tuple([gene]+list(count)))

#subject = fname[-26:][:7]
#counts = np.loadtxt(fname, skiprows=1, dtype=None, usecols=(1,))
#genes = np.loadtxt(fname, skiprows=1, dtype=np.str, usecols=(0,))

#with open('data.csv', 'wb') as csvfile:
#    writer = csv.writer(csvfile)
#    writer.writerow(('gene', subject))
#    for gene, count in zip(genes, counts):
#        writer.writerow((gene, str(count)))
