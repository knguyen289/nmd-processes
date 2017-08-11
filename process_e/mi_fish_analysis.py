import pandas as pd
import copy

from kn_tools.basic_tools import text_to_df

mi_df = text_to_df('./another_try/mutual_inf.tsv')

mi_df['MI'] = map(float,mi_df['MI'])
mi_df['Fisher'] = map(float,mi_df['Fisher'])

analy_df = copy.deepcopy(mi_df)
analy_df.insert(len(analy_df.columns),'MI_big',0)
analy_df.insert(len(analy_df.columns),'Fish_small',0)
analy_df.insert(len(analy_df.columns),'bs',0)

for index,row in analy_df.iterrows():
    if float(row.get_value('MI')) > 0.7:
        analy_df.loc[index,'MI_big'] = 1
    if float(row.get_value('Fisher')) < 0.05:
        analy_df.loc[index,'Fish_small'] = 1
    analy_df.loc[index,'bs'] = analy_df.loc[index,'MI_big'] * analy_df.loc[index,'Fish_small']

analy_df.sort_values(['Fisher','MI'],ascending=[True,False],inplace=True)
analy_df.reset_index(drop=True,inplace=True)

f = open('MI_Fish_cand.txt','w')
o_str = '\n'.join(list(set(analy_df[analy_df['bs'] == 1]['name2'])))
f.write(o_str)
f.close()

analy_df.to_csv('MI_Fish_analysis.tsv',sep='\t',index=False)

