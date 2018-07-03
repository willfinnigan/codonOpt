import codonOpt
import logging
from codonOpt.AnalysisRedesignTools import GCTools

logging.basicConfig(level=logging.INFO)



class Sequence(list):

    def __str__(self):
        s = ""
        for i in range(len(self)):
            s += self[i]

        return s



codon_table = codonOpt.CodonTable(codon_tables_dir='Data/', json_file='Escherichia_coli_K12.json', low_cuttoff=0.1)
codon_table.give_all_codons_below_cuttoff_equal_weighting()
seq_gen = codonOpt.Sequence_Generator(codon_table, check_dna_to_protein=False)

protein_seq = "MRAVVFENKERVAVKEVNAPRLQHPLDALVRVHLAGICGSDLHLYHGKIPVLPGSVLGHEFVGQVEAVGEGIQDLQPGDWVVGPFHIACGTCPYCRRHQYNLCERGGVYGYGPMFGNLQGAQAEILRVPFSNVNLRKLPPNLSPERAIFAGDILSTAYGGLIQGQLRPGDSVAVIGAGPVGLMAIEVAQVLGASKILAIDRIPERLERAASLGAIPINAEQENPVRRVRSETNDEGPDLVLEAVGGAATLSLALEMVRPGGRVSAVGVDNAPSFPFPLASGLVKDLTFRIGLANVHLYIDAVLALLASGRLQPERIVSHYLPLEEAPRGYELFDRKEALKVLLVVRGGGSGDYKDDDDK**"

def make_dna(protein, seq_gen):
    dna_seq = []

    for aa in protein:
        dna_seq.append(seq_gen.optimise(aa))

    return dna_seq


"""Types"""
from deap import creator, base, tools, algorithms


creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
creator.create("Individual", Sequence, fitness=creator.FitnessMin)

toolbox = base.Toolbox()


toolbox.register("make_dna", make_dna, protein_seq, seq_gen)


toolbox.register("individual",tools.initIterate, creator.Individual, toolbox.make_dna)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

print(str(toolbox.individual()))


def evaluate(individual):
    from codonOpt.AnalysisRedesignTools import GCTools
    no_over = GCTools.percentage_bp_in_windows_over_gc_threshold(str(individual), 100, 60)

    return(no_over,)


def myMutation(individual, sequence_gen=''):
    import random
    start = random.randint(0,len(individual))
    end = random.randint(0,len(individual) - start)

    for i in range(start,end):
        individual[i] = sequence_gen.optimise_from_dna(individual[i])

    return (individual,)

toolbox.register("mate", tools.cxTwoPoint)
toolbox.register("select", tools.selBest)
toolbox.register("evaluate", evaluate)
toolbox.register("mutate", myMutation, sequence_gen=seq_gen)


if __name__ == "__main__":
    pop = toolbox.population(n=1000)

    fit = 0.0
    for i in range(1000):
        algorithms.eaMuPlusLambda (pop, toolbox,
                                   400, 100, #parents, children
                                   .2, .4, #probabilities
                                   1) #iterations

        top = sorted(pop, key=lambda x:x.fitness.values[0])[-1]
        fit = top.fitness.values[0]
        print('---------------------')
        print('Fit = ' + str(fit))
        print('Top = ' + str(top))


