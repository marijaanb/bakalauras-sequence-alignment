from asyncio import subprocess
import random
import subprocess

N = 10
ilgis = 2000000

bash = "daryk_bazes.sh"

b = open(bash, 'w')

b.write('rm seka*.fasta.* \n')

fragmentas = 'CAGCTAAGGTAGCTACCAATATTTAGTTTTTTAGCCTTGCGACAGACCTCCTACTTAGATTGCCACGCATTGAGCTAGCGAGTCAGCGATAAGCATGACG'
kiek_mutaciju = 5
def mutink(seka, n):
    atsakymas = seka
    for i in range(n):
        atsitiktine_pozicija = random.randint(0, len(seka) -1)
        raide = random.sample("ACGT", 1)[0]
        atsakymas = atsakymas[:atsitiktine_pozicija] + str(raide) + atsakymas[atsitiktine_pozicija+1:]
    return atsakymas

for i in range(N):
    s = ">seka" + str(i) + '\n'
    for j in range(ilgis):
        s = s + random.sample("ACGT", 1)[0]
    atsitiktine_pozicija = random.randint(100, len(s) - 100)
    mut_frgmentas = mutink(fragmentas, kiek_mutaciju)
    s = s[:atsitiktine_pozicija] + mut_frgmentas + s[atsitiktine_pozicija:]
    s = s+'\n'
    failas = 'seka' + str(i) + '.fasta'
    print(failas)
    f = open(failas, 'w')
    f.write(s)
    
    f.close()
    b.write(f'/Users/marijaanb/Desktop/blast/ncbi-blast-2.12.0+/bin/makeblastdb \
            -in {failas} -dbtype nucl \n')
b.close()    

    #failas = '/Users/marijaanb/Desktop/blast/ncbi-blast-2.12.0+/bin/random_db' + '/seka' + str(i) + '.fasta'
db = subprocess.run(["bash", "daryk_bazes.sh"])
    #db = subprocess.run(["head", "seka1.fasta"])