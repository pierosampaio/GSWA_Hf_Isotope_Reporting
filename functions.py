################################################################################
######################### Input processing functions ###########################
################################################################################

import numpy as np
import pandas as pd
import os
from scipy.optimize import brentq
from scipy.stats import norm
import re
from definitions import *


def filter_by_regex(strings, pattern):
  """
  Returns a list of strings from `strings` that match the given regex `pattern`.
  
  Args:
      strings (list): List of strings to search.
      pattern (str): Regular expression pattern.
  
  Returns:
      list: Strings that match the pattern.
  """
  if not isinstance(strings, list) or not all(isinstance(s, str) for s in strings):
      raise ValueError("Input must be a list of strings.")
  if not isinstance(pattern, str):
      raise ValueError("Pattern must be a string.")

  try:
      regex = re.compile(pattern)
  except re.error as e:
      raise ValueError(f"Invalid regular expression: {e}")
  
  matches = [s for s in strings if regex.search(s)]

  if len(matches) == 1: return matches[0]
  elif not matches: return None
  else: print(matches)


################################################################################
######################### Input processing functions ###########################
################################################################################

def filter_by_regex(strings, pattern):
    """
    Returns a list of strings from `strings` that match the given regex `pattern`.
    
    Args:
        strings (list): List of strings to search.
        pattern (str): Regular expression pattern.
    
    Returns:
        list: Strings that match the pattern.
    """
    if not isinstance(strings, list) or not all(isinstance(s, str) for s in strings):
        raise ValueError("Input must be a list of strings.")
    if not isinstance(pattern, str):
        raise ValueError("Pattern must be a string.")

    try:
        regex = re.compile(pattern)
    except re.error as e:
        raise ValueError(f"Invalid regular expression: {e}")
    
    matches = [s for s in strings if regex.search(s)]

    if len(matches) == 1: return matches[0]
    elif not matches: return None
    else: print(matches)


def LuHf_process(LuHf_file, Sample_list = None):
  
  """
  Ensure the input spreadsheets have the columns in the following order:
  1: Sample/spot
  2: Duration
  3: 176Hf/177Hf
  4: 176Hf/177Hf 2SE
  5: 176Lu/177Hf
  6: 176Lu/177Hf 2SE
  7: 176Yb/177Hf
  8: 176Yb/177Hf 2SE
  9: 178Hf/177Hf
  10: 178Hf/177Hf 2SE
  11: Total Hf (V)
  12: TK comment (optional) - will be created if not present
  """

  
  
  # Read all sheets into a dictionary of DataFrames
  all_sheets = pd.read_excel(LuHf_file, sheet_name=None)
  
  # Concatenate them into one DataFrame
  df = pd.concat(all_sheets.values(), ignore_index=True)
  
  # Get rid of skipped rows
  df = df.dropna(how="all",axis=0)

  df = df.set_index("Sample/spot")

  pattern = r"(?i)^duration(\s*\(s\))?$"
  duration_col = filter_by_regex(list(df.columns),pattern)

  if duration_col == None:
      duration_col = "Duration"
      df.insert(loc=0, column=duration_col, value=np.zeros(len(df)))

  # Ensure the existence of TK comment column
  if "TK comment" not in df.columns:
      placeholdercol = np.empty(len(df))
      df["TK comment"] = placeholdercol.fill(np.nan)

  df = df.loc[~df.index.str.contains("FC1|MTZ|MUN|OGC|91500")]

  
  # Create sample column
  SampleName = []
  
  for i, idx in enumerate(df.index):
      base_str = str(idx).rsplit("-",1)[0]
      pattern = r"(?i)^[a-z].*"
      if re.match(pattern, base_str):
          new_str = str(idx).rsplit("-",1)[0].split("-")[-1]
      else:
          new_str = str(idx).rsplit("-",1)[0].split("-")[0]
      
      SampleName.append(new_str)

  
  df["Sample"] = SampleName
  df = df.loc[~df.Sample.str.contains("FC1|MTZ|MUN|OGC|91500")]

      

  
  df = df.dropna(
      subset = duration_col, axis = 0
  )
  
  df = df.loc[:,~df.columns.str.startswith("Unnamed")]
  

  col_names = [
      "Duration","Hf176Hf177","Hf176Hf177_2SE","Lu176Hf177","Lu176Hf177_2SE",
      "Yb176Hf177","Yb176Hf177_2SE","Hf178Hf177","Hf178Hf177_2SE",
      "Total_Beam","TK_comment","Sample"
  ]
  

  df = df.rename(
      dict(zip(df.columns,col_names)),
      axis = 1
  )


  df = df.dropna(axis = 0, how = "all")
  #df = df.dropna(axis = 1, how = "any")
  
  

  if Sample_list:
      # Match shortened names in the analysis index to the full sample names
      # Remove extensions
      
      full_samples = [samp.split(".")[0] for samp in Sample_list]
  
      rename_dict = {}

      print(df.head())
      
    
      for df_sample in df["Sample"].unique():
          # Find all full sample names that end with this df_sample
      
          matches = [full for full in full_samples if full.endswith(df_sample)]
  
          if len(matches) == 1:
              rename_dict[df_sample] = matches[0]
          elif len(matches) > 1:
              # If multiple matches, choose the longest (most specific)
              rename_dict[df_sample] = max(matches, key=len)
          # else: no match → leave unchanged
  
      df["Sample"] = df["Sample"].replace(rename_dict)

      print("Samples in Lu-Hf dataset:")
      print("\n".join(list(df.Sample.unique())))
      print("\n")

  df["Spot"] = df.index.values
  df = df.rename_axis("id", axis = 0)
  

  df["SampleSpot"] = df.apply(
        lambda r: "-".join([
        r["Sample"],
        *r["Spot"].split("-")[1:]   
        ]),axis=1)

  df["SampleSpot"] = df["SampleSpot"].str.replace(".","-")

  

  mask = df["SampleSpot"].astype("str").str.match(r".+-\d+$") & ~df["SampleSpot"].astype("str").str.match(r".+-\d+-\d+$")

  df.loc[mask, "SampleSpot"] += "-1"

  df["SampleSpot"] = df["SampleSpot"].astype("str").str.removeprefix("GSWA_")
  df = df.rename({"Spot":"LU_HF_ANALYSIS_ID"},axis=1)
  try:
      df = df.drop("Standards")
      df = df.drop("TK_comment",axis=1)
  except:
      pass
  
  return df

def UPb_xls_process(UPb_path,UPb_file):

  header = [
      "GroupID","SpotNo","GrainSpot","238U_ppm","232Th_ppm",
      "Th232U238","f204_pct","U238Pb206","U238Pb206_1sig",
      "Pb207Pb206","Pb207Pb206_1sig","4corr_U238Pb206",
      "4corr_U238Pb206_1sig","4corr_Pb207Pb206","4corr_Pb207Pb206_1sig",
      "4corr_86_date","4corr_86_date_1sig","4corr_76_date",
      "4corr_76_date_1sig","Discordance_pct"
  ]

  Sample_name = UPb_file.split(".")[0]
  Sample_name = Sample_name.removesuffix("-combined")

  xl = pd.ExcelFile(os.path.join(UPb_path,UPb_file))
  sheet_names = xl.sheet_names
  target_sheets = ["excel_table","data_table"]
  found_sheet = next((s for s in target_sheets if s in sheet_names), None)
  skiprows = {
      "excel_table":1,
      "data_table":4
  }

  df = pd.read_excel(os.path.join(UPb_path,UPb_file), sheet_name = found_sheet, skiprows = skiprows[found_sheet])

  ## Sanitise dataframe
  df = df.iloc[:,:20]
  df = df.dropna(axis=0)


  df = df.rename(
  dict(zip(df.columns,header)),
  axis=1
  )
  
  

  if found_sheet == "data_table":
      df["Spot"] = [str(v).replace(".","-") for v in df["GrainSpot"].values]
  elif found_sheet == "excel_table":
      df["Spot"] = df["GrainSpot"].str.split("-", expand=True)[1].str.replace(".","-")

  
  df["UPB_ANALYSIS_ID"] = df.GrainSpot
  

  df["Sample"] = Sample_name
  
  df["SampleSpot"] = (df["Sample"] + "-" + df["Spot"])

  df["SampleSpot"] = df["SampleSpot"].str.replace(".","-")

  df["SpotNo"] = df["SpotNo"].astype("int").astype("string")
  
  
  
  return(df)


def UPb_txt_process(UPb_path,UPb_file):
  
  header = [
  "Sample","GroupID","SpotNo","UPB_ANALYSIS_ID","238U_ppm","232Th_ppm",
  "Th232U238","f204_pct","U238Pb206","U238Pb206_1sig",
  "Pb207Pb206","Pb207Pb206_1sig","4corr_U238Pb206",
  "4corr_U238Pb206_1sig","4corr_Pb207Pb206","4corr_Pb207Pb206_1sig",
  "4corr_86_date","4corr_86_date_1sig","4corr_76_date",
  "4corr_76_date_1sig","Discordance_pct","SampleSpot","GrainSpot"
  ]

  df = pd.read_table(
      os.path.join(
          UPb_path,UPb_file
      )
  )
  df["Geochronid"] = df["Geochronid"].astype("string")
  df["Grain spot"] = df["Grain spot"].astype("string")
  

  df["Grain spot2"] = df["Grain spot"].str.replace(".","-")
  
  if df["Geochronid"].str.startswith("GSWA").any():
      df["SampleSpot"] = (
          df["Geochronid"].str.split("_", expand = True)[1].str.split(".", expand = True)[0] +
          "-" +
          df["Grain spot2"]
      ).astype("str")
  else: 
      df["SampleSpot"] = (
          df["Geochronid"].str.split(".", expand = True)[0] + "-" +
          df["Grain spot2"]
      ).astype("str")

  ## SANITISE AND RENAME COLUMNS
  
  df = df.rename({"Grain Spot":"UPB_ANALYSIS_ID"}, axis = 1)
  
  df = df[[
      "Geochronid","Grp_ID","Spot no", "UPB_ANALYSIS_ID", "238U(ppm)",
      "232Th(ppm)", "232Th_238U", "f(%)", "238U_206Pb", "238U_206Pb_er",
      "207Pb_206Pb", "207Pb_206Pb_er", "238U_206Pb*", "238U_206Pb*_er",
      "207Pb*_206*Pb", "207Pb*_206*Pb_er", "238U_206Pb*_age",
      "238U_206Pb*_age_er", "207Pb*_206Pb*_age", "207Pb*_206Pb*_age_er",
      "Disc(%)","SampleSpot"
  ]]

  df = df.rename(
      dict(zip(df.columns,header)),
      axis = 1 
  )

  return df
    
def UPb_file_join(UPb_path):

  UPb_dfs = []
  for f in os.listdir(UPb_path):
      if os.path.splitext(f)[-1] == ".xls": 
          UPb_dfs.append(UPb_xls_process(
              UPb_path,
              f
          ))
      elif os.path.splitext(f)[-1] == ".txt":
          UPb_dfs.append(UPb_txt_process(
              UPb_path,
              f
          ))
      else:
          pass
  assert len(UPb_dfs) > 0, "No files to concatenate"
  return pd.concat(UPb_dfs)

def merge_datasets(UPb_dataset: pd.DataFrame, LuHf_dataset: pd.DataFrame, joining_key: str) -> pd.DataFrame:

  UPb_dataset[joining_key] = UPb_dataset[joining_key].astype("str")
  LuHf_dataset[joining_key] = LuHf_dataset[joining_key].astype("str")
  
  df_merged =  pd.merge(UPb_dataset, LuHf_dataset, on = joining_key)
  
  
  if len(LuHf_dataset) > len(df_merged) or len(UPb_dataset) > len(df_merged):
      Problem = list(set(LuHf_dataset[joining_key]) - set(df_merged[joining_key]))
      print("There are Hf analyses not present in the U-Pb dataset")
  
  print("Manually check the following analyses numbers:\n")
  for p in Problem:
      print(p)
      

  return df_merged





################################################################################
############################ Hf spot calculations ##############################
################################################################################



################################## Functions ###################################


def calc_initial_ratios(Hf176, Lu176, Age):
  return Hf176 - Lu176*(np.exp(Lambda_Lu*Age*1e6) - 1)

def initial_ratio_uncertainty(Hf176_unc):
  return Hf176_unc

def epsilon_Hf(Hf176, Lu176, Age):
  eps = ((calc_initial_ratios(Hf176, Lu176, Age)/calc_initial_ratios(Hf176_CHUR, Lu176_CHUR, Age)) - 1) * 1e4
  return eps

def epsilon_Hf_uncertainty(Hf176_unc, Age):
  # Assuming no uncertainty on CHUR value for simplicity
  # Taking the uncertainty of the back-calculated value
  return 1e4 * (Hf176_unc/calc_initial_ratios(Hf176_CHUR, Lu176_CHUR, Age))
    
def two_stage_TDM_YJ(Hf176, Lu176, Age):
  """
  Calculates two stage TDM model ages. Assumes first differentiation at 4500 Ma. And
  an average crustal 176Lu/177Hf of 0.012. The latter can be changed in the definitions
  file
  """
  DM_init = calc_initial_ratios(Hf176_DM,Lu176_DM,Age)
  return ((Age/1000) + ((1/(Lambda_Lu*1e9)) * np.log(1+((Hf176-DM_init))/(Lu176_Crust-Lu176_DM))))*1000


def two_stage_TDM_v2(Hf176,Lu176,Age):
  # after Janousek (2024)

  t = Age * 1e6

  numerator = Hf176 - (np.exp(Lambda_Lu*t)-1)*(Lu176-Lu176_Crust)-Hf176_DM
  denominator = Lu176_Crust - Lu176_DM
  return ((np.log(numerator/denominator+1))/Lambda_Lu)*1e-6


#####################################################################################
############################# weighted mean calculations ############################
#####################################################################################

def calculate_mswd(vals, errs, w_mean):
  n = len(vals)
  if n <=1: return 0
  weights = 1./ (errs**2)
  mswd = np.sum(weights * (vals - w_mean)**2)/(n-1)
  return mswd


def weighted_mean(values, uncertainties, method = "internal_sigma"):
  """
  Performs a weighted mean calculation with inverse variance weighting. Optionally
  performs outlier rejection.

  Parameters:
  - values: data points
  - uncertainties: 1-sigma uncertainties (SE or SD)
  - method: "chauvenet" (as in IsoplotR), "internal_sigma" (GSWA
  age calculations) or None (does not exclude outliers)
  """

  vals = np.array(values, dtype=float)
  errs = np.array(uncertainties, dtype=float)

  iteration = 0
  while True:
      n = len(vals)
      if n <= 1: break

      weights = 1.0/(errs**2)

      w_mean = np.dot(vals,weights)/np.sum(weights)

      # Calculate deviations
      deviations = np.abs(vals - w_mean)

      if method.lower() == "chauvenet":
      # spread-based: uses Z-score (deviation/population_std)
          std_dev = np.std(vals)
          if std_dev == 0: break
          z_scores = deviations/std_dev
          probs = n * (2 * (1 - norm.cdf(z_scores)))

          if np.min(probs) < 0.05:
              reject_idx = np.argmin(probs)
              do_reject = True
          else:
              do_reject = False

      elif method.lower() == "internal_sigma":
          # Uncertainty-based: deviation/individual error
          # Checks if the mean is > X sigma away from the point's own error
          sigma_ratios = deviations/errs
          max_ratio_idx = np.argmax(sigma_ratios)

          if sigma_ratios[max_ratio_idx] > 2.5:
              reject_idx = max_ratio_idx
              do_reject = True
          else:
              do_reject = False

      if do_reject:
          vals = np.delete(vals, reject_idx)
          errs = np.delete(errs, reject_idx)
          iteration += 1
      else:
          break

  final_w_mean = np.sum((1.0/(errs**2)) * vals) / np.sum(1.0/(errs**2))
  internal_err = np.sqrt(1.0/np.sum(1.0/(errs**2)))
  mswd = calculate_mswd(vals,errs,final_w_mean)
  external_err = internal_err * np.sqrt(mswd)
  ci = 1.96 * internal_err


  return {
      "mean": final_w_mean,
      "uncertainty": internal_err,
      "external_error": np.round(external_err,1),
      "ci": np.round(ci,1),
      "mswd": mswd,
      "data": vals,
      "n_rejected": iteration
  }
    


def calc_group_stats(df: pd.DataFrame, grouping_var: str, variable: str, variable_unc: str, prepend: str) -> pd.DataFrame: 
  """
  Takes a grouped dataframe and calculates group level aggregate statistics for a given
  variable, which should be the name of the column of interest and its associated uncertainty
  (variable_unc)
  """
  
  
  
  dfg = df.groupby(grouping_var)
  
  median = dfg[variable].median()
  median_1sigma = np.sqrt((dfg[variable_unc].median())**2 + (dfg[variable].std())**2) # Following Yong Jun's method
  average = dfg[variable].mean()
  std = dfg[variable].std()
  maxx = dfg[variable].max()
  min = dfg[variable].min()
  
  sample_weighted_mean = (
          dfg
          .apply(lambda g: weighted_mean(
              g[variable],
              g[variable_unc]
          )
      ))
  
  sample_weighted_mean = pd.DataFrame.from_dict(sample_weighted_mean.to_dict(), orient="index")
  wm = sample_weighted_mean["mean"]
  uncertainty = sample_weighted_mean["uncertainty"]
  
  data = np.array((
      median,median_1sigma,average,std,maxx,min,wm,uncertainty
  ))
  
  headers = [
      "MDN","1SIG","AVE","1SD","MAX","MIN","WEIGHTEDMEAN","WEIGHTEDMEAN_95_pct_CI"
  ]
  
  index = sample_weighted_mean.index
  columns = [prepend + "_" + lab for lab in headers]
  
  return pd.DataFrame(data=data.T,columns=columns,index=index)


def overwrite_igneous_age(df: pd.DataFrame) -> pd.DataFrame:
  
  dff = df.copy(deep=True)
  
  print("Available samples:")
  for sample in dff['Sample'].unique():
    print(f"\n{sample}")
  
  sample = ""
  while sample.upper() != "X":
    sample = input("Select sample to replace age. \nType 'X' to finish editing ages. Type 'Exit' to quit the function.\nSample:")
    if sample.upper() == "X": break
    if sample.title() == "Exit": return None
    
    print(f"Current age for {sample}: {dff.loc[np.logical_and(dff['Sample']==sample,dff.GroupID == 'I'),'Age_Hf_calculation'].values[0]: .1f} \
    ± {dff.loc[np.logical_and(dff['Sample']==sample,dff.GroupID == 'I'),'Age_Hf_calculation_unc'].values[0]: .1f}")
    new_age_unc_raw = input("Type new age and uncertainty in the format XXXX/XX:")
    new_age, new_unc = float(new_age_unc_raw.split("/")[0]),float(new_age_unc_raw.split("/")[1])
    
    dff.loc[np.logical_and(dff["Sample"] == sample, dff.GroupID == "I"),
    "Age_Hf_calculation"] = new_age
    
    dff.loc[np.logical_and(dff["Sample"] == sample, dff.GroupID == "I"),
    "Age_Hf_calculation_unc"] = new_unc

  return dff


def create_aggregate_df(df: pd.DataFrame, grouping_var: str, vals: list, uncertainties: list, prepends: list, by_list = None) -> pd.DataFrame:

  df.insert(len(df.columns),"exclude",np.nan)
  
  
  
  if by_list is not None:
    with open(by_list, "r") as f:
      vals_to_exclude = [line.strip() for line in f if line.strip() and not line.startswith("#")]
  
  
  else:
    print("Select spots to exclude:")
    vals_to_exclude = []
    
    for val in df[grouping_var].unique():
        print(f"\n\nSpots in sample {val}:\n")
  
        available_spots = df.loc[df[grouping_var] == val, "SampleSpot"].tolist()
        print("\n".join(available_spots))
        
        print("Select spot number to exclude from calculation. Type just the spot number, without the sample.\nType 'X' to finish for this sample. Type 'Exit' to quit the function.")
        s = ""
        while s != "X":
            s = input("Spot number:")
            if s.upper() == "X": break
            if s.title() == "Exit": return None
            if (val.split(".")[0] + "-" + s) not in map(str,available_spots):
                print("Invalid spot!")
                continue
            vals_to_exclude.append(val.split(".")[0] + "-" + s)
        corr = input("Confirm selection? Press 'Y' to confirm or 'N' to select spots to remove from list.")
        if corr.upper() == "N":
          print("Select spot undo, press X to finish:")
          s = ""
          while s.upper() != "X":
            s = input("Spot number:")
            if s.upper() == "X": break
            if s.title() == "Exit": return None
            if (val.split(".")[0] + "-" + s) not in vals_to_exclude:
              print("Valid not in exclude list.")
              continue
            vals_to_exclude.remove(val.split(".")[0] + "-" + s)
            
  print("Excluded spots:\n","\n".join(vals_to_exclude))
          

  df.loc[df["SampleSpot"].isin(vals_to_exclude),"exclude"] = 1

  df = df.loc[df["exclude"].isna()]

  dataframes = []
  for v,u,p in zip(vals, uncertainties, prepends):
      dataframes.append(calc_group_stats(df,grouping_var,v,u,p))

  master_df = pd.concat(dataframes, axis = 1)
  
  
  
  try:
    ages = df.loc[df.GroupID == "I"].groupby("Sample")["Age_Hf_calculation"].mean()
    sigs = df.loc[df.GroupID == "I"].groupby("Sample")["Age_Hf_calculation_unc"].mean()
    
    master_df["Age"] = master_df.index.map(ages)
    master_df["Age_1sig"] = master_df.index.map(sigs)
    
  except:
    pass
  
  return master_df


### U-Pb constants for calculating weighted mean ages

lbd235 = 9.8485e-10
lbd238 = 1.55125e-10
U238_U235 = 137.818

def ratio_function(t):
  return (1/U238_U235) * \
         (np.exp(lbd235*t) - 1) / \
         (np.exp(lbd238*t) - 1)

def age_from_ratio(ratio):
  def f(t):
      return ratio_function(t) - ratio
  return brentq(f, 1e6, 4.5e9)

def pb207_pb206_age_with_uncertainty(ratio, sigma_ratio):
    
  # Solve for age
  t = age_from_ratio(ratio)
  
  # Compute derivative dR/dt
  A = np.exp(lbd235*t)
  B = np.exp(lbd238*t)
  
  dRdt = (1/U238_U235) * (
      (lbd235*A*(B-1) - lbd238*B*(A-1)) /
      ((B-1)**2)
  )
  
  # Propagate uncertainty
  sigma_t = sigma_ratio / abs(dRdt)
  
  t *= 1e-6
  sigma_t *= 1e-6

	
  return t, sigma_t
