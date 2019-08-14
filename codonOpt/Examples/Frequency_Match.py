import codonOpt
import math

seq = codonOpt.DNA("ATGCGCGCCGTTGTATTTGAGAATAAAGAACGCGTCGCCGTGAAAGAAGTGAACGCACCTCGTCTGCAGCACCCGTTAGATGCTCTCGTTCGCGTCCATCTGGCGGGGATTTGTGGCTCTGATTTACACCTGTATCACGGTAAAATCCCAGTGCTGCCGGGAAGCGTGCTGGGTCACGAATTCGTCGGGCAAGTTGAAGCCGTGGGCGAAGGCATCCAGGATCTGCAGCCGGGCGATTGGGTGGTCGGTCCATTCCATATTGCGTGTGGAACGTGTCCTTATTGTCGTCGCCACCAGTACAATCTTTGTGAGCGTGGTGGCGTTTATGGTTACGGCCCGATGTTTGGCAACCTCCAAGGAGCGCAAGCAGAAATCCTCCGCGTCCCTTTTTCGAATGTGAACTTACGCAAACTGCCACCGAACCTCAGCCCCGAGCGCGCAATTTTTGCGGGCGACATTCTCAGCACAGCCTATGGCGGACTCATCCAGGGCCAGCTGCGCCCTGGTGATTCAGTAGCGGTAATTGGTGCAGGCCCAGTGGGGCTGATGGCGATCGAAGTGGCGCAGGTGTTGGGCGCGTCAAAAATCCTGGCGATTGACCGTATCCCGGAGCGTCTGGAACGCGCTGCGTCGCTGGGCGCCATCCCAATTAATGCAGAACAGGAAAACCCGGTTCGTCGCGTTCGCAGTGAAACCAACGATGAAGGCCCCGACTTGGTGCTGGAAGCTGTAGGTGGAGCGGCGACACTTTCTCTGGCGCTCGAGATGGTTCGCCCTGGTGGTCGCGTTTCTGCAGTCGGCGTAGACAACGCACCGAGCTTCCCATTTCCACTGGCCTCGGGATTGGTTAAAGATCTCACGTTTCGTATCGGACTGGCGAATGTGCATCTGTATATTGACGCAGTTCTGGCCTTACTCGCGAGTGGTCGCCTTCAGCCGGAACGCATTGTCAGTCATTATTTGCCACTGGAAGAGGCCCCGCGCGGTTACGAACTGTTCGATCGCAAAGAGGCGCTGAAGGTGCTGCTGGTGGTGCGTGGCGGTGGGTCCGGTGATTACAAGGACGACGATGATAAATGATAA")
print(seq)
print(seq.protein)

codon_table_source = codonOpt.CodonTable(codon_tables_dir='Examples/Data/', json_file='Escherichia_coli_K12.json')
seq.get_frequencies(codon_table_source)
print(seq.freq)

codon_table_uniform = codonOpt.CodonTable()
codon_table_uniform.make_uniform()
seq_gen = codonOpt.Sequence_Generator(codon_table_uniform)

new_seq = codonOpt.DNA(codonOpt.make_dna(seq.protein, seq_gen))

codon_table_destination = codonOpt.CodonTable(codon_tables_dir='Examples/Data/', json_file='Tth codon table.json')

def evaluate(test_dna, source_dna=seq, codon_table_destination=codon_table_destination):
    differences = []

    test_dna.get_frequencies(codon_table_destination)

    for i in range(len(test_dna.freq)):
        test_freq = test_dna.freq[i]
        source_freq = source_dna.freq[i]

        diff = test_freq - source_freq
        differences.append(math.sqrt(diff*diff))

    return (sum(differences),)

print(evaluate(new_seq))

"""Types"""
from deap import creator, base, tools, algorithms


creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", codonOpt.DNA, fitness=creator.FitnessMin)

toolbox = base.Toolbox()

toolbox.register("make_dna", codonOpt.make_dna, seq.protein, seq_gen)
toolbox.register("individual",tools.initIterate, creator.Individual, toolbox.make_dna)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

print(str(toolbox.individual()))

def myMutation(individual, sequence_gen=''):
    import random
    start = random.randint(0,len(individual))
    end = random.randint(0,len(individual) - start)

    for i in range(start,end):
        individual[i] = sequence_gen.optimise_from_dna(individual[i])

    return (individual,)

ind = toolbox.individual()
print(ind)
ind = myMutation(ind, sequence_gen=seq_gen)[0]
print(ind)
ind = myMutation(ind, sequence_gen=seq_gen)[0]
print(ind)
ind = myMutation(ind, sequence_gen=seq_gen)[0]
print(ind)
ind = myMutation(ind, sequence_gen=seq_gen)[0]
print(ind)
ind = myMutation(ind, sequence_gen=seq_gen)[0]
print(ind)

toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("select", tools.selBest)
toolbox.register("evaluate", evaluate)
toolbox.register("mutate", myMutation, sequence_gen=seq_gen)


if __name__ == "__main__":
    mu = 400 #number for next gen
    lamb = 100 #children in next gen
    cxpb = 0.2 #The probability that an offspring is produced by crossover.
    mutpb = 0.4 #The probability that an offspring is produced by mutation
    ngen = 1


    pop = toolbox.population(n=1000)

    fit = 0.0
    for i in range(1000):
        algorithms.eaMuPlusLambda (pop, toolbox,
                                   mu, lamb, #parents, children
                                   cxpb, mutpb, #probabilities
                                   ngen) #iterations

        top = sorted(pop, key=lambda x:x.fitness.values[0])[-1]
        fit = top.fitness.values[0]
        print('---------------------')
        print('Fit = ' + str(fit))
        print('Top = ' + str(top))
        print("Src = " + str(seq.freq))
        print("Top = " + str(top.freq))




