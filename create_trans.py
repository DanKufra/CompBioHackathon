import numpy as np
import csv

AVG_NUM_M = 6
AVG_LEN_M = 160
AVG_LEN_I = 18459.46
PROB_STAY = 0.2
PROB_LIVE = 0.8
M_BASNAT = 3
I_BASNAT = 7

def create_trans(number_m_s, number_i_s, is_end, letters_num, states):
    states_num = 1 if is_end else 0
    states_num+=number_i_s+number_m_s
    my_t = np.zeros((states_num,states_num))
    if is_end:
        last_row = np.array([0]*(states_num-1)+[1])
        my_t[-1]= last_row
    for i in range(min(M_BASNAT,number_m_s)):
        my_row = np.zeros(states_num)
        my_row[i] = PROB_STAY
        my_row[i+1] = PROB_LIVE
        my_t[i] = my_row
    for i in range(M_BASNAT,number_m_s):
        state = i+1
        my_row = np.zeros(states_num)
        intron_transit_p = max(state/(159.56-letters_num*state), state/159.56)
        my_row[number_m_s] = intron_transit_p
        p_left = 1-intron_transit_p
        my_row[i] = PROB_STAY*p_left
        my_row[i+1] = PROB_LIVE*p_left
        my_t[i] = my_row
    for i in range(number_m_s, number_m_s+min(I_BASNAT,number_i_s-1)):
        my_row = np.zeros(states_num)
        my_row[i] = PROB_STAY
        my_row[i+1] = PROB_LIVE
        my_t[i] = my_row
    for i in range(I_BASNAT+number_m_s,number_m_s+number_i_s-1):
        my_row = np.zeros(states_num)
        state = i-number_m_s+1
        exon_transit_p = max(6/(18459.46-(letters_num*state*100)),6/18459.46)
        my_row[number_m_s-1] = exon_transit_p
        p_left = 1-exon_transit_p
        my_row[i] = PROB_STAY*p_left
        my_row[i+1] = PROB_LIVE*p_left
        my_t[i] = my_row
    last_row = np.zeros(states_num)
    if is_end:
        last_row[-2] = 0.9
        last_row[-1] = 0.1
    else:
        last_row[-1] = 1.
    my_t[number_m_s+number_i_s-1] =last_row
    with open('output.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter=' ')
        tsv_writer.writerow(['index']+list(states))
        for i in range(states_num):
            tsv_writer.writerow([states[i]] + list(my_t[i]))

create_trans(10,10,0,2,["A"]*20)