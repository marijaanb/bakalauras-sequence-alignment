from itertools import zip_longest, product
from random import choice
import subprocess 
import operator
import numpy as np
import random
import pandas as pd
import matplotlib.pyplot as plt
from dotenv import load_dotenv
import os


load_dotenv()
POPULATION_PATH=os.getenv('POPULATION_PATH')
BLAST_PATH=os.getenv('BLAST_PATH')
REAL_DB_PATHS=os.getenv('REAL_DB_PATH')
ARTIFICIAL_DB_PATHS=os.getenv('ARTIFICIAL_DB_PATHS')
OUTPUT_PATH=os.getenv('OUTPUT_PATH')


class Individual:
    def __init__(self, nucleotide):
        self.nucleotide = nucleotide

    def __repr__(self) -> str:
        return f"'{self.nucleotide}'"
    
    def __len__(self):
        return len(self.nucleotide)


class Simulation:
    def __init__(self, my_db, crossover_strategy, mutation_strategy, score_reduction, sample_size=5,
                 seq_length = 10, crossover_rate = 50, mutation_rate=50,
                 max_iter=100):
        self.sample_size = sample_size
        self.nucleotides = ["G", "A", "C", "T"]
        self.seq_length = seq_length
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate
        self.population = []
        # self.exit_threshold = exit_threshold
        self.population_path = POPULATION_PATH
        self.my_db = my_db
        self.crossover_strategy = crossover_strategy
        self.mutation_strategy = mutation_strategy
        self.best_score_per_iteration = {}
        self.max_iter = max_iter
        self.score_reduction = score_reduction
    
    def sample(self, n):
        for i in range(n):
            nucleotide=""
            for count in range(self.seq_length):
                nucleotide+=choice(self.nucleotides)
            self.population.append(Individual(nucleotide))
        self.population[3] = 'CAGCTAAGGTAGCTACCAATATTTGGTTTTTTAGCCTTGCGACAGACCTCCTACTTAGACTGCCACGCATTGAGCTAGCGAGTCAGCGATTAGCATGACG'
    
    def simulate(self):
        self.sample(self.sample_size)
        
        i = 0
        while True:
            i += 1
            # write population to FASTA file for BLAST
            print("-"*10, f"Iteration {i}", "-"*10)
            with open(self.population_path, "w") as population_file:
                for p in range(len(self.population)):
                    self.population[p] = str(self.population[p]).replace("'","").replace("'","")
                    print('>' + str(p), '\n' + self.population[p], file=population_file)
            
            self.selection(i)
            if i >= self.max_iter:
                return self.best_score_per_iteration
                

    def fitness(self):
        
        my_db = self.my_db
        all_seq_score_dicts = {}
        # run blast
        for i in range(len(my_db)):
            get_seq_score = subprocess.run([f"{BLAST_PATH} -query {POPULATION_PATH} \
                                        -db {my_db[i]} \
                                        -word_size 13\
                                        -outfmt '6 qseqid bitscore'", "/dev/null"], shell = True, capture_output=True)
            # parse blast output
            decoded = get_seq_score.stdout.decode("utf-8")
            chunks = decoded.split('\n')

            # dictionary with query id as key and max score as value for each sequence that got a score
            seq_score_dict = {}
            for c in chunks:
                if c != "":
                    idx, val = c.split('\t') if '\t' in c else c.split('\n')
                    seq_score_dict[int(idx)] = max(seq_score_dict.get(int(idx), 0), float(val))
                    
            # database dict
            all_seq_score_dicts[my_db[i]] = seq_score_dict
        
        # max score from db for each sequence
        if self.score_reduction == 'sum':
            max_scores = pd.DataFrame(all_seq_score_dicts).sum(axis=1).to_dict()
        elif self.score_reduction == 'min':
            max_scores = pd.DataFrame(all_seq_score_dicts).min(axis=1).to_dict()
        sorted_max_scores = dict(sorted(max_scores.items(), key=operator.itemgetter(1),reverse=True))
        
        ind = [*sorted_max_scores]
        
        # sort population
        self.population = [self.population[i] for i in ind]

        return self.population, sorted_max_scores
        
    def selection(self, iteration):
        # get sorted population from fitness func
        sorted_population, sc = self.fitness()
        self.best_score_per_iteration[iteration] = list(sc.values())[0]
        
        print(f"Best score: {list(sc.values())[0]}, Best sequence: {sorted_population[0]}, Last sequence: {sorted_population[-1]}")
        
        # find the number of strong individuals (20% of population)
        strong_seq_num = int(len(self.population) * 0.25)
        
        # select strong individuals
        strong_sequences = sorted_population[:strong_seq_num]
        self.population = strong_sequences
        
        # crossover
        newborn = self.crossover(strong_sequences)
        
        # add newborns to population
        self.population.extend(newborn)
        
    
    def crossover(self, strong_sequences):
        newborn = []
        
        for n in range(self.sample_size - len(self.population)): 
            seq_1 = random.choice(strong_sequences)
            seq_2 = random.choice(strong_sequences)
            
            # single point crossover
            if self.crossover_strategy == 'single_point':
                cut_seq_1 = seq_1[:len(seq_1)//2]
                cut_seq_2 = seq_2[len(seq_2)//2:]
                new_nucleotide = cut_seq_1 + cut_seq_2
                new_sequence = Individual(self.mutation(new_nucleotide))
                newborn.append(new_sequence)
            
            # uniform
            elif self.crossover_strategy == 'uniform':
                new_nucleotide = ''
                for base_1, base_2 in zip_longest(seq_1, seq_2):
                    dice = random.randint(1, 100)
                    if dice >= self.crossover_rate and base_1 is not None and base_2 is not None:
                        new_nucleotide += random.choice([base_1, base_2])
                    elif base_2 is not None and base_1 is None:
                        new_nucleotide += base_2
                    else: 
                        new_nucleotide += base_1
            
                new_sequence = Individual(self.mutation(new_nucleotide))
                newborn.append(new_sequence)
            
        return newborn

    def mutation(self, sequence):
        
        dice = random.randint(1, 100)
        # prolongation
        if self.mutation_strategy == 'prolongation':
            mutated_seq = list(sequence)
            base = choice('ATGC')
            if self.mutation_rate < dice:
                mutated_seq = [base] + mutated_seq
            else:
                mutated_seq = mutated_seq + [base]
            return ''.join(mutated_seq)
            
        
        # inversion mutation
        elif self.mutation_strategy == 'inversion':
            if self.mutation_rate < dice:
                substr_len = random.randint(2, len(sequence))
                idx = random.randrange(0, len(sequence) - substr_len + 1)
                mutated_seq = sequence[:idx] + sequence[idx : (idx+substr_len)][::-1] + sequence[(idx+substr_len):]
                return mutated_seq
        
        # swap mutation
        elif self.mutation_strategy == 'swap':
            mutated_seq = list(sequence)
            if self.mutation_rate < dice:
                for i in range(15):
                    k = random.sample(range(1, len(mutated_seq)-1), 2)
                    mutated_seq[k[0]], mutated_seq[k[1]] = mutated_seq[k[1]], mutated_seq[k[0]]
                mutated_seq = ''.join(mutated_seq)
                return mutated_seq
            
        # scramble mutation
        elif self.mutation_strategy == 'scramble':
            if self.mutation_rate < dice:
                substr_len = random.randint(2, len(sequence))
                idx = random.randrange(0, len(sequence) - substr_len + 1)
                sh = list(sequence[idx : (idx+substr_len)])
                random.shuffle(sh)
                sh = ''.join(sh)
                mutated_seq = sequence[:idx] + sh + sequence[(idx+substr_len):]
                return mutated_seq
            
        return sequence
        
if __name__ == '__main__':
    # ------------------- Define params -------------------
    mutation_rates = [0, 20, 50]
    databases = [
        ARTIFICIAL_DB_PATHS.split(":"),
        REAL_DB_PATHS.split(":")
        ]
    db_mapping = {
        'artificial': databases[0],
        'real': databases[1]
    }
    db_names = ['artificial', 'real']
    crossover_strategies = ['single_point', 'uniform']
    mutation_strategies = ['prolongation', 'inversion', 'swap', 'scramble']
    max_iter = 50
    score_reductions = ['sum', 'min']
    # ------------------- Create results DataFrame -------------------
    columns = list(range(1, max_iter+1))+['mutation_rate', 'db', 'crossover', 'mutation', 'reduction']
    results_dataframe = pd.DataFrame(columns=columns)
    # ------------------- Run simulations and save results to a file -------------------
    products = product(mutation_rates, db_names, crossover_strategies, mutation_strategies, score_reductions)
    for n_experiment, (mutation_rate, db, crossover, mutation, reduction) in enumerate(products):
        params = {
            n_experiment: 
                {
                "mutation_rate": mutation_rate,
                "db": db,
                "crossover": crossover,
                "mutation": mutation,
                "reduction": reduction
                }
            }
        print(params)
        simulation = Simulation(my_db=db_mapping[db],
                                crossover_strategy=crossover, 
                                mutation_strategy=mutation, 
                                score_reduction=reduction,
                                sample_size=50,
                                seq_length=100,
                                crossover_rate=50, 
                                mutation_rate=mutation_rate, 
                                max_iter=max_iter)
        scores = simulation.simulate()
        scores = {n_experiment: scores}
        
        df_entry = pd.concat([pd.DataFrame(scores).T, pd.DataFrame(params).T], axis=1)
        results_dataframe = pd.concat([results_dataframe, df_entry])
        
        results_dataframe.to_csv(OUTPUT_PATH)
