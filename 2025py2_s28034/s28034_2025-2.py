from Bio import Entrez
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
        
def search_and_fetch(email, api_key, tax_id, min_seq_len, max_seq_len, start=0, max_records=10):
    try:
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'
        search_term = f"txid{tax_id}[Organism] AND {min_seq_len}:{max_seq_len}[Sequence Length]"
        search_handle = Entrez.esearch(
            db='nucleotide',
            term=search_term,
            retstart=start,
            retmax=max_records,
            usehistory='n'
        )
        search_results = Entrez.read(search_handle)
        id_list = search_results['IdList']
        print(id_list)
        if not id_list:
            print("No records found matching criteria.")
            return []
        fetch_handle = Entrez.efetch(
            db='nucleotide',
            id=",".join(id_list),
            retmode='xml'
        )
        records_xml = Entrez.read(fetch_handle)
        return records_xml
    except Exception as e:
        print(f'Exception occurred:\n{e}')
        return []

def export_to_csv(records, filename="genbank_records.txt"):
    data = []
    for record in records:
        accession = record.get('GBSeq_primary-accession', 'N/A')
        seq_len = record.get('GBSeq_length', 'N/A')
        seq_desc = record.get('GBSeq_definition', 'N/A')
        
        data.append({
            "Accession Number": accession,
            "Sequence Length": seq_len,
            "Description": seq_desc
        })
    df = pd.DataFrame(data)
    df.to_csv('genbank_data.csv')

def plot_sequence_lengths(csv_file):
    df = pd.read_csv(csv_file)
    df['Sequence Length'] = pd.to_numeric(df['Sequence Length'], errors='coerce')
    df = df.dropna(subset=['Sequence Length'])
    df = df.sort_values(by='Sequence Length', ascending=False)
    plt.figure(figsize=(12, 6))
    sns.lineplot(data=df, x='Accession Number', y='Sequence Length', marker='o')
    plt.xticks(rotation=45, ha='right')
    plt.title('Sequence Lengths by Accession Number (Descending)')
    plt.tight_layout()
    plt.show()

def test():
    samples_xml = search_and_fetch(
        email='s28034@pjwstk.edu.pl', 
        api_key='5d333e4eefd7225207577a2acf13d1498008',
        tax_id=9615,
        min_seq_len=20,
        max_seq_len=100,
    )
    export_to_csv(samples_xml)

if __name__ == "__main__":
    # main()
    # test()
    plot_sequence_lengths('genbank_data.csv')