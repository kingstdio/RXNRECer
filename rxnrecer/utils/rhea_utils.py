'''
Author: Zhenkun Shi
Date: 2023-06-28 15:14:34
LastEditors: Zhenkun Shi
LastEditTime: 2023-07-03 12:25:23
FilePath: /preaction/tools/rheatool.py
Description: 

Copyright (c) 2023 by tibd, All Rights Reserved. 
'''


import re




def get_chebi_smiles(chebi_accession, chebi_cmp_df):
    match_item = chebi_cmp_df[chebi_cmp_df.chebi_accession==chebi_accession]
    if len(match_item)>0:
        res = match_item.smiles.values[0]
        if res == None:
            return '*'
        else:
            return res
    else:
        return '*'


def format_equation_chebi(rheaid,equation, chebiid, chebi_cmp_df):
    if '=' not in equation:
        print(equation)
        return ''
    equation = equation.replace(' carrying a second [4Fe-4S](2+) cluster','').replace('trans,octa-cis','trans-octa-cis')

    equation_array= equation.split(' = ') #Split reactants and products
    substrates_array = equation_array[0].split(' + ') # substrates
    products_array = equation_array[1].split(' + ') # products
    chebiid_array = chebiid.split(';') #Split chebi_id

    eq_len = len(substrates_array) + len(products_array) # Number of compounds in equation
    chebi_len = len(chebiid_array)

    # If chemical equation and chebiid counts match, expand once
    if eq_len == chebi_len:
        eq_chebi = ' + '.join(chebiid_array[0:len(substrates_array)]) + ' = ' + ' + '.join(chebiid_array[len(substrates_array):])

    else: # Has circulating metabolites, compound counts don't match
        # regex_str = r'\d |\(in\)|\(out\)|n[+-]?\d* '
        # regex_str = r'[\d ]|\(in\)|\(out\)|n|\(n[+-]?\d?\)|L-leucyl-'
        regex_str = r'\d |\(in\)|\(out\)|n|\(n[+-]?\d?\)|L-leucyl-'

        substrates_array = [re.sub(regex_str, '', item) for item in substrates_array] # Replace coefficients
        products_array = [re.sub(regex_str, '', item) for item in products_array] # Replace coefficients

        sub_products_array = []
        pos_counter = 0
        try:
            if (len(substrates_array)==len(products_array)) and (len(products_array)==len(chebiid_array)):
                sub_products_array = chebiid_array
            else:
                for i in range(len(products_array)):        # Process one compound at a time
                    if products_array[i] in substrates_array:
                        sub_products_array.append(chebiid_array[substrates_array.index(products_array[i])])
                    else:
                        sub_products_array.append(chebiid_array[(len(substrates_array))+pos_counter])
                        pos_counter+=1
            eq_chebi = ' + '.join(chebiid_array[0:len(substrates_array)]) + ' = ' + ' + '.join(sub_products_array)
        except Exception as e:
            print(rheaid)
            print(equation)
            # print(substrates_array)
            # print(products_array)
            print(e)
            eq_chebi=''
            str_equ_smiles = ''
    
    try:
        chebi_substrates_array = eq_chebi.split(' = ')[0].split(' + ')
        chebi_products_array = eq_chebi.split(' = ')[1].split(' + ')
    
        str_substrates_smile = '.'.join( [get_chebi_smiles(chebi_accession=item, chebi_cmp_df=chebi_cmp_df) for item in chebi_substrates_array]) # Reaction SMILES string
        str_products_smile = '.'.join( [get_chebi_smiles(chebi_accession=item, chebi_cmp_df=chebi_cmp_df) for item in chebi_products_array]) # Reaction SMILES string
        str_equ_smiles = f'{str_substrates_smile}>>{str_products_smile}'
    except Exception as e:
        print(rheaid)
        print(equation)
        print(eq_chebi)
        print(e)
        str_equ_smiles = ''

    return eq_chebi, str_equ_smiles