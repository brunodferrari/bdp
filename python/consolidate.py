# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 18:54:24 2021

@author: Bruno Ferrari
"""

_path="C:/Users/Bruno Ferrari/Documents/Bruno/2019/2s/MC/artigos revisão/Artigos Mes/GD/"

import pandas as pd


df_vns = pd.DataFrame()
df_ts = pd.DataFrame()
df_gs_vns = pd.DataFrame()
df_gs = pd.DataFrame()


df_vns = pd.concat([pd.read_excel(_path+"bdp/dbdp_instances/meta_vns.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_vns_1.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_vns_2.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_vns_3.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_vns_4.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_vns_5.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_vns_6.xlsx", sheet_name="results"), 
                        ], axis=0).set_index("Instance")

df_ts = pd.concat([pd.read_excel(_path+"bdp/dbdp_instances/meta_tabu.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_tabu_1.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_tabu_2.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_tabu_3.xlsx", sheet_name="results"),
                        ], axis=0).set_index("Instance")

df_gs =  pd.concat([pd.read_excel(_path+"bdp/dbdp_instances/meta_grasp_1_new_cross.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_grasp_2_new_cross.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_grasp_3_new_cross.xlsx", sheet_name="results"),
                        ], axis=0).set_index("Instance")

df_gs_vns =  pd.concat([pd.read_excel(_path+"bdp/dbdp_instances/meta_grasp2_1.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/meta_grasp2_2.xlsx", sheet_name="results"),
                        ], axis=0).set_index("Instance")

df_results= pd.concat([pd.read_excel(_path+"bdp/dbdp_instances/metafeat.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/metafeat2.xlsx", sheet_name="results"),
                        pd.read_excel(_path+"bdp/dbdp_instances/metafeat3.xlsx", sheet_name="results")
                        ], axis=0).set_index("Instance")

df_final = df_results.join([df_vns[['Crossing', 'Time']].add_suffix('_vns'), 
                 df_ts[['Crossing', 'Time']].add_suffix('_ts'), 
                 df_gs[['Crossing', 'Time']].add_suffix('_gs'), 
                 df_gs_vns[['Crossing', 'Time']].add_suffix('_gs_vns')])


df_final.columns = ['V', 
                    'E',
                    'V_diff',
                    'deg_mean',
                    'deg_std',
                    'deg_median',
                    'density',
                    'V1',
                    'V2',
                    'deg_min']+list(df_final.columns[10:])

with pd.ExcelWriter(_path+"bdp/dbdp_instances/final_results.xlsx") as writer:
    df_final.to_excel(writer, sheet_name="results", index=True)