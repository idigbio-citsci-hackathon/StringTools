import re
from collections import deque
from sys import argv

class NameSplitter:
    
    def __init__(self):
        pass
        
    def split(self, input):
        '''Given a list of collectors as a string, try to parse it into an array of individual names.'''
        
        # Preliminary clean-up:
    
        input = re.sub(r'^collected by\b:?', '', input, flags=re.I)
        input = re.sub(r'\b(collectors?)|(collrs?\.?)$', '', input, flags=re.I)
        input = re.sub(r'(?<![A-Z])([A-Z]),', r'\1.', input) # There are lots of commas after initials.
        input = re.sub(r'\.+', '.', input)
        
        # Before tokenization, look at the very end and see if there are trailing initials.
        # This would suggest that the inputs are transposed: If they are, try to switch them back. 
        # I haven't seen that many instances where input are inverted in this data,
        # so the match is pretty strict, since it could really confuse things is misapplied.
        
        strict_init_re = r'(?:[A-Z]\.\s?){1,2}'
        strict_last_re = r'\w{2,}'
        between_re = r',\s+'
        
        if re.search(strict_last_re + between_re + strict_init_re + '$', input):
            input = re.sub(
                r'\b(' + strict_last_re + ')' + between_re + '(' + strict_init_re + ')(?=,|$)',
                r'\2 \1', input)
        
        # Assuming first names precede last, split into tokens:

        tokens = deque(re.split(
            r'''
                (?:  \s+ \W? (?:with|and) \W? \s+ ) |
                (?:  \bw/ ) |
                (?:  [,;+&()/] )                        # not confident that splitting on parens is best.
            ''', input, flags=re.X|re.I))
        
        names = []
        while tokens:
            token = tokens.popleft()
            if names and re.match(
                    r'''
                        ^\s* (
                            jr\.? |
                            sr\.? |
                            ph\.?d\.?
                        ) \s* ,? \s* $
                    ''', token, flags=re.X|re.I):
                names[len(names) - 1] += ', ' + token
            elif token:
                names.append(token)
                
        names = [re.sub(r'\s+', ' ', name.strip(' ')) for name in names]
        
        # Distribute last names across first names:
        #   Scan from right-to-left;
        #     if there's a last name, store it
        #     if it's just a first name, and we have a last name, append it. 
        
        full_names = deque()
        first_or_init_re = r'''(?:
            (?: \w{2,} \.? )                    # name, possibly followed by period. (too fragile?)
            | (?: (?: [A-Za-z] \.? \s? ){1,2} ) # initials w/o periods.
        )''' # This might break either with short last names, or JRR Tolkein.
        last_name = ''
        
        while names:
            current = names.pop()
            match = re.search(
                first_or_init_re + r'''
                    \s+
                    ([\w-]{2,}) # last name
                    $ # no trailing punctuation
                ''', current, re.X)
            if match:
                last_name = match.group(1)
            elif last_name and re.match('^' + first_or_init_re + '$', current, re.X):
                current += ' ' + last_name
            if current:
                full_names.appendleft(current)
        
        return list(full_names)
        
        
if __name__ == '__main__':
    splitter = NameSplitter()
    try:
        print(splitter.split(argv[1]))
    except IndexError:
        print(NameSplitter.split.__doc__)