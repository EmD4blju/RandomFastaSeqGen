from Bio import Entrez 
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
        
def search_and_fetch(email, api_key, tax_id, min_seq_len, max_seq_len, start=0, max_records=10):
    Entrez.email, Entrez.api_key, Entrez.tool = email, api_key, 'BioScriptEx10'
    ids = Entrez.read(Entrez.esearch(db='nucleotide', term=f"txid{tax_id}[Organism] AND {min_seq_len}:{max_seq_len}[Sequence Length]", retstart=start, retmax=max_records, usehistory='n'))['IdList']
    if ids:
        records = Entrez.read(Entrez.efetch(db='nucleotide', id=",".join(ids), retmode='xml'))
        return len(ids), records
    return 0, []

if __name__ == "__main__":
    obtained_records_num, samples_xml = search_and_fetch(email=input('Input your email:'), api_key=input('Input api key:'), tax_id=input('Input taxonomy id:'), min_seq_len=input('Minimal sequence length:'),max_seq_len=input('Maximal sequence length:'))
    print(f'Obtained {obtained_records_num} records.')
    file_path = Path('saved_samples_gb', f'{input('Input filename:')}.csv')
    pd.DataFrame([{"Accession Number": record.get('GBSeq_primary-accession', 'N/A'),"Sequence Length": record.get('GBSeq_length', 'N/A'),"Description": record.get('GBSeq_definition', 'N/A')} for record in samples_xml]).to_csv(file_path) 
    sns.lineplot(data= pd.read_csv(file_path).sort_values(by='Sequence Length', ascending=False), x='Accession Number', y='Sequence Length', marker='o')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()