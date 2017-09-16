
import subprocess

# Warn if Vienna not found.
from shutil import which
if not which("RNAfold"):
    from sys import stderr
    print("WARNING: Vienna RNA is not installed.  Please install.", file=stderr)
del(which)

class RNA_fold_output():

    def __init__(self, sequence, program="RNAfold"):

        rnafold = subprocess.Popen([program, "--noPS"],
                                   stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   # Universal Newlines effectively allows string IO.
                                   universal_newlines=True)

        rna_fold_output, rna_fold_error = rnafold.communicate(str(sequence))

        if rna_fold_error:
            raise NameError('Vienna RNA Error')

        output_lines = rna_fold_output.strip().splitlines()

        self.structure = output_lines[1].split(None, 1)[0].strip()
        self.energy = float(output_lines[1].rsplit("(", 1)[1].strip("()").strip())
        self.sequence = output_lines[0]









if __name__ == "__main__":
    seq1 = "ACCTCCTTTAGAGAGAGAGAGAGAGGGGGGGGGCCCCC"
    seq2 = "TGGAGGAAA"
    seq = seq1 + "&" + seq2

    fold = RNA_fold_output(seq1, program='RNAfold')

    print(fold.structure)
    print(fold.energy)