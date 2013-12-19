## DEFINE TEST STRING

x = ["The quick brown fox", "The quick Brown fox", "Tha quick brown fox jumps"]

import string
import nltk
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import subprocess
import re
from Bio import AlignIO
from Bio.Align import AlignInfo


## APPROACH ONE: TOKEN ALIGNMENT
def token_consensus(x, wdir):
    ''' Arguments:
        x       -- List of strings

        Returns:
        Consensus string
    '''
    # Concatenate strings and tokenize 
    regexptoken = nltk.RegexpTokenizer(pattern = '\s+', gaps = True)                           
    y = string.join(x)
    tokens = regexptoken.tokenize(y)
    # tokens = nltk.word_tokenize(y) #word_tokenize assumes you are a sentence, so any ultimate periods are tokenize separately

    # Find unique tokens
    unique_tokens = []
    for token in tokens:
        if token not in unique_tokens:
            unique_tokens.append(token)
            
    # Create a dictionary of tokens
    unique_token_id = dict()
    index = 0
    for unique_token in unique_tokens:
        unique_token_id[unique_token] = string.ascii_letters[index]
        index += 1

    # Convert strings into token IDs
    entry_token_strings = []
    for entry in x:
        entry_tokens = regexptoken.tokenize(entry) #nltk.word_tokenize(entry)
        entry_token_string = ""
        for token in entry_tokens:
            #print token
            entry_token_string += unique_token_id[token]
            #print entry_token_string
        entry_token_strings.append(entry_token_string)

    # Convert strings into SeqRecord objects
    entry_token_strings = [Seq(token) for token in entry_token_strings]
    temp = [SeqRecord(entry_token_string, id = str(acc)) for acc, entry_token_string in enumerate(entry_token_strings)]

    # Write a fasta file
    temp_file = os.path.join(wdir, "temp.fasta")
    SeqIO.write(temp, temp_file, "fasta")

    # Align using alignment algorithm MAFFT
    mafft = '/usr/local/bin/mafft'
    res = subprocess.check_output([mafft, '--text', temp_file])
    out_file = os.path.join(wdir, "temp_align.fasta")
    with file(out_file, 'w') as f:
        f.write(res)

    alignres = AlignIO.read(out_file, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignres)

    # Determine consensus
    consensus = summary_align.dumb_consensus(threshold = 0.5, consensus_alpha = "alphabet", require_multiple = 1, ambiguous = "")

    # Reinterpret consensus
    inv_unique_token_id = {v:k for k, v in unique_token_id.items()}
    consensus = [inv_unique_token_id[token] for token in consensus]
    consensus = string.join(consensus, sep = " ")
    
    return str(consensus)

## APPROACH 2: CHARACTER ALIGNMENT
def character_consensus(x, wdir):
    ''' x   -- list of strings

        Returns: Single string
    '''

    # Convert string into SeqRecord objects
    y = [re.sub(" ", "_", string) for string in x]
    temp = [Seq(string, "alphabet") for string in y]
    temp = [SeqRecord(string, id = str(acc)) for acc, string in enumerate(temp)]
    temp_file = os.path.join(wdir, "temp.fasta")
    SeqIO.write(temp, temp_file, "fasta")

    # Subprocessing MAFFT
    mafft = '/usr/local/bin/mafft'
    res = subprocess.check_output([mafft, '--text', temp_file])

    # Export results into working dir
    out_file = os.path.join(wdir, "temp_align.fasta")
    with file(out_file, 'w') as f:
        f.write(res)

    # Import alignment
    alignres = AlignIO.read(out_file, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignres)

    # Determine consensus
    consensus = summary_align.dumb_consensus(threshold = 0.5, require_multiple = 1, ambiguous = "")

    return re.sub("_", " ", str(consensus))


## TESTING NFN TRANSCRIPTIONS AGAINST GOLD DATASET (i.e. transcribed in verbatim)
y = ["12 mi. W. Oakland", "12 mi. West Oakland", "12 mi W. Oakland", "12 miles W Oakland"]
test_dir = "/Users/junyinglim/Desktop"

asd1 = character_consensus(y, wdir = test_dir)
asd2 = token_consensus(y, test_dir)

import pandas as pd
from collections import defaultdict
import itertools

nfn_data = pd.read_csv(os.path.join(test_dir,"Calbug_NfN.csv"))
nfn_data = nfn_data.fillna("") # Converts all NaNs into empty strings

gold_data = pd.read_csv(os.path.join(test_dir,"Calbug_Gold.csv"))
gold_data = gold_data.fillna("") # Converts all NaNs into empty strings

# Generate dictionary of accessions paired with list of field entries
entry_list = zip(list(nfn_data["filename"]), list(nfn_data["Collector"]))
entry_id = defaultdict(list)
for k,v in entry_list:
    entry_id[k].append(v)


# Find consensus in NfN data
entry_results = defaultdict(list)
for k,v in entry_id.iteritems():
    print "Reconciling transcriptions for", k

    # If entries are identical, then entry is consensus
    if len(set(v)) == 1:
        entry_results[k].append(v[0])

    # If entries are not identical, use consensus
    else:
        entry_results[k].append(character_consensus(v, test_dir))

# Convert results into dataframe OMG I'm starting to love list comprehensions
collector = [str(collector[0]) for collector in entry_results.values()] # Necessary to index 0 and default dict values are lists
filename = [str(filename) for filename in entry_results.keys()]
results = pd.DataFrame({"filename":filename, "Collector":collector})

# Comparison with gold data set
compare_results = pd.merge(gold_data, results, on = "filename", suffixes = ("_gold", "_consensus"))


x = compare_results["Collector_consensus"] == compare_results["Collector_gold"]
compare_results.to_csv(os.path.join(test_dir, "prelim.csv"))


## APPROACH 3: SEQUENCE MATCHER APPROACH (might be faster?)
# Sequence align token IDs

# Search for first name until last name found
# 
# Remove tokens that are JUST punctuations

import nltk
from difflib import SequenceMatcher
import json

humantext = "Blah blah"
ocrtext = "8lah blah"

ocrsub = {}
human_tokens = nltk.word_tokenize(humantext.lower())
ocr_tokens = nltk.word_tokenize(ocrtext.lower())

s = SequenceMatcher(None,human_tokens,ocr_tokens)
ocrsub["result"] = json.dumps({
   "matches": s.get_matching_blocks(),
                   "tokens": ocr_tokens,
   "score": {
       "ratio": s.ratio()
   }
})
ocrsub["score"] = 100*s.ratio()
