## DEFINE TEST STRING

x = ["The quick brown fox", "The quick Brown fox", "Tha quick brown fox jumps"]

## APPROACH ONE: TOKEN ALIGNMENT
import string
import nltk
import os

# Generate tokens 
y = string.join(x)
tokens = nltk.word_tokenize(y)

# Find unique tokens
unique_tokens = []
for token in tokens:
    if token not in unique_tokens:
        unique_tokens.append(token)
        
# Create a dictionary of tokens
unique_token_id = dict()
index = 0
for unique_token in unique_tokens:
    unique_token_id[unique_token] = string.ascii_uppercase[index]
    index += 1

# Convert strings into token IDs
entry_token_strings = []
for entry in x:
    entry_tokens = nltk.word_tokenize(entry)
    entry_token_string = ""
    for token in entry_tokens:
        #print token
        entry_token_string += unique_token_id[token]
        #print entry_token_string
    entry_token_strings.append(entry_token_string)

# Convert strings into SeqRecord objects
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

test_dir = "/Users/junyinglim/Desktop"
entry_token_strings = [Seq(token) for token in entry_token_strings]


test = [SeqRecord(entry_token_string, id = str(acc)) for acc, entry_token_string in enumerate(entry_token_strings)]

# Write a fasta file
from Bio import SeqIO
test_file = os.path.join(test_dir, "test.fasta")
SeqIO.write(test, test_file, "fasta")

# Align using alignment algorithm MAFFT
mafft = '/usr/local/bin/mafft'
import subprocess
res = subprocess.check_output([mafft, '--anysymbol', test_file])

out_file = os.path.join(test_dir, "test_align.fasta")
with file(out_file, 'w') as f:
    f.write(res)

from Bio import AlignIO
alignres = AlignIO.read(out_file, "fasta")

from Bio.Align import AlignInfo
summary_align = AlignInfo.SummaryInfo(alignres)

# Determine consensus
consensus = summary_align.dumb_consensus(threshold = 0.5, consensus_alpha = "alphabet", require_multiple = 1, ambiguous = "")


## APPROACH 2: CHARACTER ALIGNMENT
import re
from Bio import AlignIO
from Bio.Align import AlignInfo

def character_consensus(x, wdir):
    ''' x   -- list of strings

        Returns: Single string
    '''

    # Convert string into SeqRecord objects
    y = [re.sub(" ", "_", string) for string in x]
    temp = [Seq(string, "alphabet") for string in y]
    temp = [SeqRecord(string, id = str(acc)) for acc, string in enumerate(temp)]
    temp_file = os.path.join(wdir, "temp.fasta")
    SeqIO.write(test, temp_file, "fasta")

    # Subprocessing MAFFT
    res = subprocess.check_output([mafft, '--text', test_file])

    # Export results into working dir
    out_file = os.path.join(test_dir, "temp_align.fasta")
    with file(out_file, 'w') as f:
        f.write(res)

    # Import alignment
    alignres = AlignIO.read(out_file, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignres)

    # Determine consensus
    consensus = summary_align.dumb_consensus(threshold = 0.5, require_multiple = 1, ambiguous = "")

    return [alignres, consensus, re.sub("_", " ", str(consensus))]

y = ["12 mi. W. Oakland", "12 mi. West Oakland", "12 mi W. Oakland", "12 miles W Oakland"]
test = character_consensus(y, test_dir)


## APPROACH 3: SEQUENCE MATCHER APPROACH
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
