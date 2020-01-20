import pandas as pd
import numpy as np
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.HMM.MarkovModel import MarkovModelBuilder, HiddenMarkovModel
from Bio import Alphabet
from Bio.Seq import MutableSeq

def build_mm(emission_df, transition_df):
    state_alphabet = Alphabet.Alphabet()
    state_alphabet.letters = ['F', 'M', 'L', 'I', 'E']
    emission_alphabet = Alphabet.Alphabet()
    emission_alphabet.letters = ['A', 'C', 'G', 'T', '^']
    m = MarkovModelBuilder(state_alphabet, emission_alphabet)
    for state in state_alphabet.letters:
        for state2 in state_alphabet.letters:
            prob = transition_df.iloc[state][state2]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('emission_matrix', help='emission matrix tsv')
    parser.add_argument('transition_matrix', help='transition_matrix tsv')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    state_alphabet = Alphabet.Alphabet()
    state_alphabet.letters = ['F', 'M', 'L', 'I', 'E']
    emission_alphabet = Alphabet.Alphabet()
    emission_alphabet.letters = ['A', 'C', 'G', 'T', '^']
    m = MarkovModelBuilder(state_alphabet, emission_alphabet)
    for state in state_alphabet.letters:
        for state2 in state_alphabet.letters:
            if args.emiss
    m.set_initial_probabilities({'F': 1.0})
    m.allow_transition('F', 'I', 0.01)
    m.allow_transition('F', 'F', 0.99)
    m.allow_transition('I', 'M', 0.01)
    m.allow_transition('I', 'I', 0.9)
    m.allow_transition('I', 'L', 0.09)
    m.allow_transition('M', 'M', 0.99)
    m.allow_transition('M', 'I', 0.01)
    m.allow_transition('L', 'L', 0.99)
    m.allow_transition('L', 'E', 0.01)
    m.allow_transition('E', 'E', 1.0)
    m.set_emission_score('F', 'A', 0.25)
    m.set_emission_score('F', 'C', 0.25)
    m.set_emission_score('F', 'G', 0.25)
    m.set_emission_score('F', 'T', 0.25)

    m.set_emission_score('M', 'A', 0.25)
    m.set_emission_score('M', 'C', 0.25)
    m.set_emission_score('M', 'G', 0.25)
    m.set_emission_score('M', 'T', 0.25)

    m.set_emission_score('L', 'A', 0.25)
    m.set_emission_score('L', 'C', 0.25)
    m.set_emission_score('L', 'G', 0.25)
    m.set_emission_score('L', 'T', 0.25)

    m.set_emission_score('I', 'A', 0.25)
    m.set_emission_score('I', 'C', 0.25)
    m.set_emission_score('I', 'G', 0.25)
    m.set_emission_score('I', 'T', 0.25)

    m.set_emission_score('E', '^', 1.0)
    mm = m.get_markov_model()

    cur_state = "F"
    roll_seq = MutableSeq("", state_alphabet)
    # generate the sequence
    for roll in range(100000):
        letter = np.random.choice(['A', 'C', 'G', 'T'], 1)[0]

    # add on a new roll to the sequence
    roll_seq.append(letter)
seq = roll_seq.toseq()

predicted_states, prob = mm.viterbi(seq, state_alphabet)
print('here')