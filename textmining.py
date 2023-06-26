from Bio import Entrez, Medline
import pandas as pd
import requests
from bs4 import BeautifulSoup
import time

def search_pubchem(metabolite):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{metabolite}/cids/JSON'
    response = requests.get(url)
    data = response.json()
    if 'IdentifierList' in data and 'CID' in data['IdentifierList']:
        return data['IdentifierList']['CID']
    return []

def search_metabolite_functions(metabolites):
    Entrez.email = 'your.email@example.com'
    results = []
    for metabolite in metabolites:
        # Search for metabolite functions in PubChem
        pubchem_results = search_pubchem(metabolite)
        if pubchem_results:
            results.append(pubchem_results)
        else:
            # If not found in PubChem, search in PubMed for freely available full-text papers
            term = f'{metabolite} AND "free full text"[Filter]'
            handle = Entrez.esearch(db='pubmed', sort='relevance', retmax='40', term=term)
            result = Entrez.read(handle)
            results.append(result["IdList"])
    return results

def search_google_scholar(metabolite):
    query = f'{metabolite} metabolite function'
    url = f'https://scholar.google.com/scholar?q={query}'
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'html.parser')
    result_links = soup.find_all('h3', class_='gs_rt')
    results = []
    for link in result_links:
        article_title = link.text.strip()
        article_url = link.a['href']
        article_id = extract_article_id(article_url)
        results.append(article_id)
    return results

def extract_article_id(url):
    # Extract article ID from the URL
    article_id = url.split('scholar?')[1]
    return article_id

def fetch_metabolite_details(id_lists):
    Entrez.email = 'your.email@example.com'
    results = []
    for id_list in id_lists:
        ids = ','.join(id_list)
        handle = Entrez.efetch(db='pubmed', rettype='medline', retmode='text', id=ids)
        records = Medline.parse(handle)
        records = list(records)
        functions = []
        for record in records:
            functions.append(record.get('AB', 'N/A'))
        results.append(functions)
    return results

def process_metabolites(excel_file, output_file_prefix, batch_size=5):
    # Read the Excel file
    df = pd.read_excel(excel_file)

    # Extract the metabolite names from the 'Metabolite' column
    metabolites = df['Metabolite'].tolist()

    total_metabolites = len(metabolites)
    total_batches = (total_metabolites + batch_size - 1) // batch_size  # Round up division

    for batch in range(total_batches):
        start_idx = batch * batch_size
        end_idx = (batch + 1) * batch_size
        batch_metabolites = metabolites[start_idx:end_idx]

        # Search for metabolite functions
        id_lists = search_metabolite_functions(batch_metabolites)

        # Fetch details of the metabolite functions
        results = fetch_metabolite_details(id_lists)

        # Flatten the nested list of functions
        flattened_results = [function for sublist in results for function in sublist]

        # Update the 'Function' column with the function results
        functions_dict = {metabolite: function for metabolite, function in zip(batch_metabolites, flattened_results)}
        df.loc[start_idx:end_idx-1, 'Function'] = df.loc[start_idx:end_idx-1, 'Metabolite'].map(functions_dict).fillna('N/A')

        # Save the updated DataFrame to a new Excel file
        output_file = f"{output_file_prefix}_{batch + 1}.xlsx"
        df.loc[start_idx:end_idx-1].to_excel(output_file, index=False)
        
        # Wait for 1 second before processing the next batch
        time.sleep(1)

def main():
    # Input file and output file prefix
    excel_file = "F:/NISER internship/Text-minning/unique_metabolites.xlsx"
    output_file_prefix = "F:/NISER internship/Text-minning/m"

    # Process metabolites in batches
    process_metabolites(excel_file, output_file_prefix, batch_size=100)


if __name__ == '__main__':
    main()
