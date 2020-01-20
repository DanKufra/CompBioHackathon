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


def calc_statistics(pred, lbl, pos_vals=None, neg_vals=None):
    # get overall accuracy
    overall_accuracy = np.sum(pred == lbl) / len(pred)

    # TPR
    if pos_vals  is not None:
        pred_is_pos = np.zeros(len(pred))
        for val in pos_vals:
            pred_is_pos[pred == val] = 1.0
        lbl_is_pos = np.zeros(len(lbl))
        for val in pos_vals:
            lbl_is_pos[lbl == val] = 1.0
        TPR = np.sum((pred_is_pos == lbl_is_pos) & (lbl_is_pos == 1.0)) / np.sum(lbl_is_pos)

    # TNR
    if neg_vals  is not None:
        pred_is_neg = np.zeros(len(pred))
        for val in neg_vals:
            pred_is_neg[pred == val] = 1.0
        lbl_is_neg = np.zeros(len(lbl))
        for val in neg_vals:
            lbl_is_neg[lbl == val] = 1.0
        TNR = np.sum((pred_is_neg == lbl_is_neg) & (lbl_is_neg == 1.0)) / np.sum(lbl_is_neg)

    print("ACC: %f TPR: %f TNR: %f" % (overall_accuracy, TPR, TNR))
    # get confusion matrix of preds vs labels (Multiclass)


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
            lbl_str += l
        lbl_str += 'E'
        print(lbl_str)
        intron_only = ''
        for i in range(len(lbl_str)):
            intron_only += 'I'

        intron_only = np.array(list(intron_only))
        lbl_arr = np.array(list(lbl_str))
        pred_arr = np.array(list(predicted_states))
        calc_statistics(pred_arr, lbl_arr, pos_vals=['F', 'M', 'L'], neg_vals=['I'])
