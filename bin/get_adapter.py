#!/usr/bin/env python
import sys
import gzip

def get_adapter(fastq,db):
    """A function to identify the most abundant adapter sequence in a fastq file."""
    
    # read fastq file
    with gzip.open(fastq,'rt') as f:
        file=f.read()
        query=file.split('@')
        query=['@'+seq for seq in query][1:]
        
        reads=[]
        for i,s in enumerate(query):        
          seq=s.split('\n')[1]
          seq_flattened_query = ''.join([line.strip() for line in seq])#join lines
          reads.append(seq_flattened_query)

    # read adapter sequence db
    with open(db) as d:
        database=d.read()
        database=database.split('>')
        database=['>'+seq for seq in database[1:]]#join

        adapters=[]
        for x in database:
            database_seq=x.split('\n')[1:]
            database_seq = ''.join([x.strip() for x in database_seq])#join lines
            adapters.append(database_seq)

    # match adapter sequence with fastq file
    match={}
    for i,s in enumerate(adapters):
        for j,t in enumerate(reads):
            if s in t:
                if s not in match:
                    match[s] = 1
                else:
                    match[s] +=1
                    
    # report max counts for adapter sequence
    if len(match) == 0:
        print('no adapters found... exiting')
        # sys.exit(0)
    else:
        print(f'{max(match, key=match.get)}')
    
if __name__== '__main__':
    fastq=sys.argv[1]
    db=sys.argv[2]
    get_adapter(fastq,db)
