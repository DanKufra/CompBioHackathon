import pandas as pd
import numpy as np
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.HMM.MarkovModel import MarkovModelBuilder, HiddenMarkovModel
from Bio import Alphabet
from Bio.Seq import MutableSeq
from seq_generator import getSeq


def build_mm(emission_df, transition_df, start_state):
    state_alphabet = Alphabet.Alphabet()
    state_alphabet.letters = transition_df.columns
    emission_alphabet = Alphabet.Alphabet()
    emission_alphabet.letters = emission_df.columns
    m = MarkovModelBuilder(state_alphabet, emission_alphabet)
    m.set_initial_probabilities({start_state: 1.0})
    for state in state_alphabet.letters:
        for state2 in state_alphabet.letters:
            prob = transition_df.loc[state][state2]
            m.allow_transition(state, state2, prob)
    for state in state_alphabet.letters:
        for letter in emission_alphabet.letters:
                prob = emission_df.loc[state][letter]
                m.set_emission_score(state, letter, prob)

    mm = m.get_markov_model()
    return mm


def run_viterbi(markov_model, seq, emission_array, state_array):
    state_alphabet = Alphabet.Alphabet()
    state_alphabet.letters = state_array
    emission_alphabet = Alphabet.Alphabet()
    emission_alphabet.letters = emission_array

    mutable_seq = MutableSeq("", emission_alphabet)
    # generate the sequence
    for letter in seq:
        # if letter == 'N':
        #     continue
        # add on a new roll to the sequence
        mutable_seq.append(letter[0].upper())
    mutable_seq = mutable_seq.toseq()

    predicted_states, prob = markov_model.viterbi(mutable_seq, state_alphabet)
    return predicted_states, prob


if __name__ == '__main__':
    # args = parse_args()
    emission_df = pd.read_csv('./initial_emissions_matrix.tsv', sep=' ', index_col=0)
    transition_df = pd.read_csv('./initial_transitions_matrix.tsv', sep=' ', index_col=0)

    mm = build_mm(emission_df, transition_df, 'F')

    gene_gen = getSeq('./hg19.genes.NR.chr19.exonCount2_29.bed', './hg19.2bit')
    for i, (seq, lbl) in enumerate(gene_gen):
        seq = seq.tolist()
        seq.append('^')
        predicted_states, prob  = run_viterbi(mm, seq, emission_df.columns, transition_df.columns)
        print(predicted_states)
        lbl_str = ''
        for l, c in zip(lbl, seq):
            # if c == 'N':
            #     continue
            lbl_str += l
        lbl_str += 'E'
        print(lbl_str)
        intron_only = ''
        for i in range(len(lbl_str)):
            intron_only += 'I'
        ACC = np.sum(np.array(list(lbl_str)) == np.array(list(predicted_states))) / len(lbl_str)
        ACC_intron = np.sum(np.array(list(lbl_str)) == np.array(list(intron_only))) / len(lbl_str)
        if ACC > 0.8:
            print("WOW")
        print("ACC is %f, ACC_intron is %f" % (ACC, ACC_intron))

        # break
