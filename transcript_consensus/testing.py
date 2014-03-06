# TESTING SUITE FOR CONSENSUS TOOLS
# Author: Junying Lim
# Date: 13th Jan 2014

from consensus_tools import *

# DEPENDENCIES
import pandas as pd
from collections import defaultdict
import itertools

# DIRECTORIES
test_dir = os.path.join(os.getenv("HOME"), "Desktop")

# TEST EXAMPLE
test_example = ["12 mi. W. Oakland, Cal", "12 mi West Oakland", "12 mi. W. Oakland", "12 miles W Oakland"]
token_example = token_consensus(test_example, test_dir)
character_example = character_consensus(test_example, wdir = test_dir)
'''
# EVALUATING NFN TRANSCRIPTIONS AGAINST GOLD DATASET (i.e. transcribed in verbatim)
nfn_data = pd.read_csv(os.path.join(test_dir,"Calbug_NfN.csv"))
nfn_data = nfn_data.fillna("") # Converts all NaNs into empty strings

gold_data = pd.read_csv(os.path.join(test_dir,"Calbug_Gold.csv"))
gold_data = gold_data.fillna("") # Converts all NaNs into empty strings

# Finding consensus
character_coll = variant_consensus(accession = "filename", field = "Collector", method = "character", data = nfn_data, wdir = test_dir)
character_loc = variant_consensus(accession = "filename", field = "Locality", method = "character", data = nfn_data, wdir = test_dir)

token_coll = variant_consensus(accession = "filename", field = "Collector", method = "token", data = nfn_data, wdir = test_dir)
token_loc = variant_consensus(accession = "filename", field = "Locality", method = "token", data = nfn_data, wdir = test_dir)

# Compiling results
character_results = pd.merge(gold_data, character_coll, on = "filename", suffixes = ("_gold", "_consensus"))
character_results = pd.merge(character_results, character_loc, on = "filename", suffixes = ("_gold", "_consensus"))

token_results = pd.merge(gold_data, token_coll, on = "filename", suffixes = ("_gold", "_consensus"))
token_results = pd.merge(token_results, token_loc, on = "filename", suffixes = ("_gold", "_consensus"))
'''
# Some summary statistics
x1 = character_results["Collector_consensus"] == character_results["Collector_gold"]
x2 = character_results["Locality_consensus"] == character_results["Locality_gold"]

y1 = token_results["Collector_consensus"] == token_results["Collector_gold"]
y2 = token_results["Locality_consensus"] == token_results["Locality_gold"]

# Export results
token_results.to_csv(os.path.join(test_dir, "token_consensus.csv"))
character_results.to_csv(os.path.join(test_dir, "character_consensus.csv"))

