## CONSENSUS TOOLS
# Description: A set of tools for reconciling multiple variants of a string of text using sequence alignment and consensus algorithms
# Author: Junying Lim
# Date: 13th Jan 2014

# Notes:
# Developed during the CITScribe hackathon organized by iDigBio, Gainsville, Florida (15 Dec - 20 Dec 2013)
# For more information on using MAFFT for string matching http://mafft.cbrc.jp/alignment/software/textcomparison.html


## DEPENDENCIES
import os # Path tools
import string # String tools
import nltk # For tokenizing
import subprocess # For subprocessing MAFFT
import re # Regular expressions
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Align import AlignInfo

mafft = "/usr/local/bin/mafft"


def token_consensus(x, wdir):
    ''' Implements a sequence alignment and consensus finding on tokenized strings

        Arguments:
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
    y = [re.sub("\.", "%", string) for string in y]
    y = [re.sub("\(", "[", string) for string in y]
    y = [re.sub("\)", "]", string) for string in y]
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
    consensus = summary_align.dumb_consensus(threshold = 0.5, require_multiple = 1, consensus_alpha = None, ambiguous = "")

    consensus = re.sub("_", " ", str(consensus))
    consensus = re.sub("%", ".", str(consensus))
    consensus = re.sub("\[", "(", str(consensus))
    consensus = re.sub("\]", ")", str(consensus))
    consensus = consensus.strip()
    return consensus

def variant_consensus(accession, field, data, method, wdir):
    
    # Generate dictionary of accessions paired with list of field entries
    entry_list = zip(list(data[accession]), list(data[field]))
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
            if method == "character":
                entry_results[k].append(character_consensus(v, wdir))
            elif method == "token":
                entry_results[k].append(token_consensus(v, wdir))


    # Convert results into dataframe [OMG I'm starting to love list comprehensions]
    est = [str(est[0]) for est in entry_results.values()] # Necessary to index 0 and default dict values are lists
    acc = [str(acc) for acc in entry_results.keys()]
    results = pd.DataFrame({str(accession):acc, str(field):est})

    # Export
    return results


# Pulled from wrapper for debugging
'''
entry_list = zip(list(nfn_data["filename"]), list(nfn_data["Locality"]))
entry_id = defaultdict(list)
for k,v in entry_list:
    entry_id[k].append(v)
entry_id["UMMZI212852 Sympetrum madidum.jpg"]
'''
## TO DO ## Convert spaces into tokens
