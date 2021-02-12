# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 10:51:07 2021

@author: REMI
"""

import sqlite3
import pandas as pd


with sqlite3.connect("../results/RNANet.db") as connection:
    chain_list = pd.read_sql("""SELECT chain_id, structure_id, chain_name
                             FROM chain JOIN structure 
                             ON chain.structure_id = structure.pdb_id
                             WHERE resolution < 4.0 
                             ORDER BY structure_id ASC
                             LIMIT 10;""",
                        con=connection)

    req = """SELECT index_chain, old_nt_resnum, nt_position, nt_name, nt_code, nt_align_code, 
                is_A, is_C, is_G, is_U, is_other, freq_A, freq_C, freq_G, freq_U, freq_other, dbn,
                paired, nb_interact, pair_type_LW, pair_type_DSSR, alpha, beta, gamma, delta, epsilon, zeta, epsilon_zeta,
                chi, bb_type, glyco_bond, form, ssZp, Dp, eta, theta, eta_prime, theta_prime, eta_base, theta_base,
                v0, v1, v2, v3, v4, amplitude, phase_angle, puckering 
                FROM 
                (SELECT chain_id, rfam_acc from chain WHERE chain_id = {})
                NATURAL JOIN re_mapping
                NATURAL JOIN nucleotide
                NATURAL JOIN align_column;"""
        
    for chain in chain_list.iterrows():
        df = pd.read_sql(req.format(chain[1].chain_id), connection)
        filename = chain[1].structure_id + '-' + chain[1].chain_name + '.csv'
        df.to_csv(filename, float_format="%.2f", index=False)            