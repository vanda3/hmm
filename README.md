# Hidden Markov Model
The model is based on Escherichia Coli.

## Model
The model separates genes from non-genes. Inside each gene we have:
- The start codon s1-s6 (ATA, ATC, ATT, ATG, GTG, CTG and TTG)
- The middle codon c1-c3 (that can be anything and can be repeated)
- The stop codon e1-e5 (TAA, TAG and TGA)
Inside each non-gene we have all four bases that transition between them. States sG, fG, sNG and fNG are “null” states, they don’t emit.

The file build.txt contains the main structure for this model.

## Algorithms
- Viterbi outputs the path for a given sequence following the most likely path of the model. However, it gets stuck in a local maximum as the protein always enters a non-gene and has a low chance of transitioning.
- Forward outputs the probability matrix.
I didn’t implement the training version of Viterbi and Forward.

## Run Code
Compile: python3 HMM.py

Then insert the sequence you want viterbi/forward to evaluate

Note: it starts by building the model, so it can take a few seconds.


