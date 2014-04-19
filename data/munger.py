import csv
import os  
import numpy as np

def get_path(name, base=os.path.dirname(__file__)):
    return os.path.join(base, name)

kirc = get_path('kirc')
kirp = get_path('kirp')
counts = {}
genes = {}

for f in os.listdir(kirc):
     ff = get_path(f, kirc)
     if os.path.isfile(ff):
        subject = ff[-26:][:7]  ## this appears to be a unique identifier, but need to verify
        print subject
        counts[subject] = np.loadtxt(ff, skiprows=1, dtype=np.str, usecols=(1,))
        genes[subject] = np.loadtxt(ff, skiprows=1, dtype=np.str, usecols=(0,))

# files = [get_path(f, kirc) for f in os.listdir(kirc)]
# header = np.array(['gene'] + [f[-26:][:7] for f in files]).reshape(74,1)
# a = [np.loadtxt(f, skiprows=1, dtype=np.str, usecols=(0,1)) for f in files]
# d = [c[:,1] for c in a if all(c[:,0] == a[0][:,0])]
# d = [a[0][:,0]] + d
# x = np.hstack((header,np.asarray(d))).T
# np.savetxt(x)

## dt = [('gene','|S25'), (subject,'|S25')]
## np.loadtxt(ff, skiprows=1, dtype=dt, usecols=(0,1))
## count.dtype.names[1]

with open('kirc.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    subjects = counts.keys()
    writer.writerow(tuple(['gene']+subjects))
    data  = np.array(counts.values())
    for gene, count in zip(genes.values()[0], data.T):
        writer.writerow(tuple([gene]+list(count)))

counts = {}
genes = {}
for f in os.listdir(kirp):
     ff = get_path(f, kirp)
     if os.path.isfile(ff):
        subject = ff[-26:][:7]  ## this appears to be a unique identifier, but need to verify
        print subject
        counts[subject] = np.loadtxt(ff, skiprows=1, dtype=np.str, usecols=(1,))
        genes[subject] = np.loadtxt(ff, skiprows=1, dtype=np.str, usecols=(0,))

with open('kirp.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    subjects = counts.keys()
    writer.writerow(tuple(['gene']+subjects))
    data  = np.array(counts.values())
    for gene, count in zip(genes.values()[0], data.T):
        writer.writerow(tuple([gene]+list(count)))
