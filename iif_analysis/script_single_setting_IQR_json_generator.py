# import sys
# sys.path.append('../../')
import numpy as np
import pandas as pd
import json
import copy
from plot_helper import coloring_legend, df_col_replace
from constant import REPORTDAYS, HEADER_NAME, COLUMNS_TO_DROP, FIRST_ROW_AFTER_BURNIN

def single_setting_IQR_json_generator(fpath_pattern_list, outfile_dir, outfile_stgy_tag, threshold):
  
  def T01_IQR_reporter_oneset_mostdangtriple(dflist, pattern, threshold):
    all_100_T01s = [] # in days
    for onerun in dflist:
      combined_geno_freq = onerun.filter(regex=pattern, axis=1).sum(axis=1).values # len=361
      T01_this_run = float('inf')
      for idx,val in enumerate(combined_geno_freq):
        if val > threshold:
          T01_this_run = REPORTDAYS[idx]
          break
      all_100_T01s.append(T01_this_run)
    assert(len(all_100_T01s)==100)
    return np.quantile(all_100_T01s, [0.25, 0.5, 0.75])

  def T01_IQR_reporter_oneset_mostdangdouble(dflist_arg, drug, threshold):
    option = 1
    most_dang_double_tag = '2-2' if drug == 'DHA-PPQ' else '2-4'
    all_100_T01s = [] # in days
    dflist = copy.deepcopy(dflist_arg)
    # rename all 100 df's by `drug` and sum-up columns
    for i in range(len(dflist)):
      dflist[i] = df_col_replace(dflist[i], drug, option)
      combined_geno_freq = dflist[i][most_dang_double_tag].values # len=361
      T01_this_run = float('inf')
      for idx,val in enumerate(combined_geno_freq):
        if val > threshold:
          T01_this_run = REPORTDAYS[idx]
          break
      all_100_T01s.append(T01_this_run)
    assert(len(all_100_T01s)==100)
    return np.quantile(all_100_T01s, [0.25, 0.5, 0.75])

  # Main Driver Code
  set3_fpath, set4_fpath, set7_fpath, set8_fpath, set11_fpath, set12_fpath = fpath_pattern_list
  # all rows, all sets
  iqr_median = {}
  iqr_25p = {}
  iqr_75p = {}

  dflist_set3 = []
  dflist_set4 = []
  dflist_set7 = []
  dflist_set8 = []
  dflist_set11 = []
  dflist_set12 = []

  for i in range(1,101):
    dflist_set3.append(
      pd.read_csv(set3_fpath%i, index_col=False, names=HEADER_NAME, sep='\t').drop(columns=COLUMNS_TO_DROP)
    )
    dflist_set4.append(
      pd.read_csv(set4_fpath%i, index_col=False, names=HEADER_NAME, sep='\t').drop(columns=COLUMNS_TO_DROP)
    )
    dflist_set7.append(
      pd.read_csv(set7_fpath%i, index_col=False, names=HEADER_NAME, sep='\t').drop(columns=COLUMNS_TO_DROP)
    )
    dflist_set8.append(
      pd.read_csv(set8_fpath%i, index_col=False, names=HEADER_NAME, sep='\t').drop(columns=COLUMNS_TO_DROP)
    )
    dflist_set11.append(
      pd.read_csv(set11_fpath%i, index_col=False, names=HEADER_NAME, sep='\t').drop(columns=COLUMNS_TO_DROP)
    )
    dflist_set12.append(
      pd.read_csv(set12_fpath%i, index_col=False, names=HEADER_NAME, sep='\t').drop(columns=COLUMNS_TO_DROP)
    )

  # initialize with row1
  # set3
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set3, 'TYY..Y2.', threshold)
  assert(len(temp)==3) # 25p, median, and 75p values
  iqr_median['row1'] = [temp[1]]
  iqr_25p['row1'] = [temp[0]]
  iqr_75p['row1'] = [temp[2]]
  # set4
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set4, 'TYY..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row1'].append(temp[1])
  iqr_25p['row1'].append(temp[0])
  iqr_75p['row1'].append(temp[2])
  # set7
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set7, 'TYY..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row1'].append(temp[1])
  iqr_25p['row1'].append(temp[0])
  iqr_75p['row1'].append(temp[2])
  # set8
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set8, 'TYY..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row1'].append(temp[1])
  iqr_25p['row1'].append(temp[0])
  iqr_75p['row1'].append(temp[2])
  # set11
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set11, 'TYY..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row1'].append(temp[1])
  iqr_25p['row1'].append(temp[0])
  iqr_75p['row1'].append(temp[2])
  # set12
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set12, 'TYY..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row1'].append(temp[1])
  iqr_25p['row1'].append(temp[0])
  iqr_75p['row1'].append(temp[2])

  # row2
  # set3
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set3, 'KNF..Y2.', threshold)
  assert(len(temp)==3) # 25p, median, and 75p values
  iqr_median['row2'] = [temp[1]]
  iqr_25p['row2'] = [temp[0]]
  iqr_75p['row2'] = [temp[2]]
  # set4
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set4, 'KNF..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row2'].append(temp[1])
  iqr_25p['row2'].append(temp[0])
  iqr_75p['row2'].append(temp[2])
  # set7
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set7, 'KNF..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row2'].append(temp[1])
  iqr_25p['row2'].append(temp[0])
  iqr_75p['row2'].append(temp[2])
  # set8
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set8, 'KNF..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row2'].append(temp[1])
  iqr_25p['row2'].append(temp[0])
  iqr_75p['row2'].append(temp[2])
  # set11
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set11, 'KNF..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row2'].append(temp[1])
  iqr_25p['row2'].append(temp[0])
  iqr_75p['row2'].append(temp[2])
  # set12
  temp = T01_IQR_reporter_oneset_mostdangtriple(dflist_set12, 'KNF..Y2.', threshold)
  assert(len(temp)==3)
  iqr_median['row2'].append(temp[1])
  iqr_25p['row2'].append(temp[0])
  iqr_75p['row2'].append(temp[2])

  # row3
  # set3
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set3, 'DHA-PPQ', threshold)
  assert(len(temp)==3) # 25p, median, and 75p values
  iqr_median['row3'] = [temp[1]]
  iqr_25p['row3'] = [temp[0]]
  iqr_75p['row3'] = [temp[2]]
  # set4
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set4, 'DHA-PPQ', threshold)
  assert(len(temp)==3)
  iqr_median['row3'].append(temp[1])
  iqr_25p['row3'].append(temp[0])
  iqr_75p['row3'].append(temp[2])
  # set7
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set7, 'DHA-PPQ', threshold)
  assert(len(temp)==3)
  iqr_median['row3'].append(temp[1])
  iqr_25p['row3'].append(temp[0])
  iqr_75p['row3'].append(temp[2])
  # set8
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set8, 'DHA-PPQ', threshold)
  assert(len(temp)==3)
  iqr_median['row3'].append(temp[1])
  iqr_25p['row3'].append(temp[0])
  iqr_75p['row3'].append(temp[2])
  # set11
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set11, 'DHA-PPQ', threshold)
  assert(len(temp)==3)
  iqr_median['row3'].append(temp[1])
  iqr_25p['row3'].append(temp[0])
  iqr_75p['row3'].append(temp[2])
  # set12
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set12, 'DHA-PPQ', threshold)
  assert(len(temp)==3)
  iqr_median['row3'].append(temp[1])
  iqr_25p['row3'].append(temp[0])
  iqr_75p['row3'].append(temp[2])

  # row4
  # set3
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set3, 'ASAQ', threshold)
  assert(len(temp)==3) # 25p, median, and 75p values
  iqr_median['row4'] = [temp[1]]
  iqr_25p['row4'] = [temp[0]]
  iqr_75p['row4'] = [temp[2]]
  # set4
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set4, 'ASAQ', threshold)
  assert(len(temp)==3)
  iqr_median['row4'].append(temp[1])
  iqr_25p['row4'].append(temp[0])
  iqr_75p['row4'].append(temp[2])
  # set7
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set7, 'ASAQ', threshold)
  assert(len(temp)==3)
  iqr_median['row4'].append(temp[1])
  iqr_25p['row4'].append(temp[0])
  iqr_75p['row4'].append(temp[2])
  # set8
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set8, 'ASAQ', threshold)
  assert(len(temp)==3)
  iqr_median['row4'].append(temp[1])
  iqr_25p['row4'].append(temp[0])
  iqr_75p['row4'].append(temp[2])
  # set11
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set11, 'ASAQ', threshold)
  assert(len(temp)==3)
  iqr_median['row4'].append(temp[1])
  iqr_25p['row4'].append(temp[0])
  iqr_75p['row4'].append(temp[2])
  # set12
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set12, 'ASAQ', threshold)
  assert(len(temp)==3)
  iqr_median['row4'].append(temp[1])
  iqr_25p['row4'].append(temp[0])
  iqr_75p['row4'].append(temp[2])

  # row5
  # set3
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set3, 'AL', threshold)
  assert(len(temp)==3) # 25p, median, and 75p values
  iqr_median['row5'] = [temp[1]]
  iqr_25p['row5'] = [temp[0]]
  iqr_75p['row5'] = [temp[2]]
  # set4
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set4, 'AL', threshold)
  assert(len(temp)==3)
  iqr_median['row5'].append(temp[1])
  iqr_25p['row5'].append(temp[0])
  iqr_75p['row5'].append(temp[2])
  # set7
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set7, 'AL', threshold)
  assert(len(temp)==3)
  iqr_median['row5'].append(temp[1])
  iqr_25p['row5'].append(temp[0])
  iqr_75p['row5'].append(temp[2])
  # set8
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set8, 'AL', threshold)
  assert(len(temp)==3)
  iqr_median['row5'].append(temp[1])
  iqr_25p['row5'].append(temp[0])
  iqr_75p['row5'].append(temp[2])
  # set11
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set11, 'AL', threshold)
  assert(len(temp)==3)
  iqr_median['row5'].append(temp[1])
  iqr_25p['row5'].append(temp[0])
  iqr_75p['row5'].append(temp[2])
  # set12
  temp = T01_IQR_reporter_oneset_mostdangdouble(dflist_set12, 'AL', threshold)
  assert(len(temp)==3)
  iqr_median['row5'].append(temp[1])
  iqr_25p['row5'].append(temp[0])
  iqr_75p['row5'].append(temp[2])

  # if directory exist check happening
  # in main script notebook file
  with open(outfile_dir+outfile_stgy_tag+'_median.json', 'w') as outfile:
    json.dump(iqr_median, outfile)
  with open(outfile_dir+outfile_stgy_tag+'_25p.json', 'w') as outfile:
    json.dump(iqr_25p, outfile)
  with open(outfile_dir+outfile_stgy_tag+'_75p.json', 'w') as outfile:
    json.dump(iqr_75p, outfile)
