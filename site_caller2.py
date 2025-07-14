import pysam
import sys

def getcall(samfile, the_chr, the_pos, ref_N, alt_N):
    cov = 0
    ac = {ref_N: 0, alt_N: 0, 'NA': 0, 'del': 0}
    iter = samfile.pileup(region=f"{the_chr}:{the_pos}-{the_pos}", stepper='nofilter')
    
    for pileupcolumn in iter:
        if pileupcolumn.reference_pos == int(the_pos) - 1:
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del:
                    ac['del'] += 1
                else:
                    obsN = pileupread.alignment.query_sequence[pileupread.query_position]
                    cov += 1
                    ac[obsN] = ac.get(obsN, 0) + 1
    
    if cov == 0 or ac["del"] > 1 or (ac["NA"] / cov) > 0.5:
        return "./."
    elif (ac[ref_N] / cov) >= 0.9:
        return "0/0"
    elif (ac[alt_N] / cov) >= 0.9:
        return "1/1"
    elif (ac[ref_N] / cov) > 0.1 and (ac[alt_N] / cov) > 0.1:
        return "0/1"
    else:
        return "./."

infile = sys.argv[1]
insites = sys.argv[2]

samfile = pysam.AlignmentFile(infile, "rb")
with open(insites, "r") as sites:
    results = []
    for line in sites:
        inputsite = line.split()
        result = getcall(samfile, inputsite[0], inputsite[1], inputsite[2], inputsite[3])
        results.append(result)

samfile.close()

print("\n".join(results))
