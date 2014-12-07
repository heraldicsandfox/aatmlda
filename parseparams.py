f = open("kickparams", "r")

uthreshs = {}
vthreshs = {}
kthreshs = {}

for line in f:
    if line.startswith("uthresh"):
        toks = line.split()
        uthresh = float(toks[2][:-1])
        vthresh = float(toks[5][:-1])
        kthresh = float(toks[8])
    elif line.startswith("Log"):
        toks = line.split()
        total_ll = float(toks[3][:-1])
        per_ll = float(toks[6])
        if uthresh != 0.1:
            uthreshs[uthresh] = (total_ll, per_ll)
        elif vthresh != 0.2:
            vthreshs[vthresh] = (total_ll, per_ll)
        elif kthresh != 0.05:
            kthreshs[kthresh] = (total_ll, per_ll)
        else:
            uthreshs[uthresh] = (total_ll, per_ll)
            vthreshs[vthresh] = (total_ll, per_ll)
            kthreshs[kthresh] = (total_ll, per_ll)

print uthreshs
print vthreshs
print kthreshs

