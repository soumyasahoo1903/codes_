import pandas as pd
from bs4 import BeautifulSoup
import requests
from Bio import Entrez
import re

# Set up NLTK and download necessary resources
import nltk
nltk.download('punkt')

# Read the Excel file
excel_file = 'F:/NISER internship/Text-minning/metabolite.xlsx'
df = pd.read_excel(excel_file)

# Extract the metabolite names from the 'Metabolite name' column
metabolites = df['Metabolites'].tolist()

# Lists to store functions, links, and hits
functions = []
links = []
hits = []

# Keywords related to the gut microbiome or gut microorganisms
gut_keywords = ['barrier integrity', 'epithelial integrity', 'intestinal stem cells', 'epithelial cell proliferation', 'wnt signalling pathway', 'cell differentiation', 'intestinal epithelial cells']

# Initialize sheet for metabolites without gut keyword mentions
no_keyword_df = pd.DataFrame(columns=['Metabolites', 'Function'])

# Iterate through each metabolite
for metabolite in metabolites:
    function = ''
    link = ''
    hit_count = 0

    # Search Metabolomics Workbench for metabolite functions
    mw_keywords = metabolite + ' ' + ' OR '.join(gut_keywords)
    mw_search_url = f"https://www.metabolomicsworkbench.org/search/?query={mw_keywords}"
    mw_response = requests.get(mw_search_url)
    if mw_response.status_code == 200:
        soup = BeautifulSoup(mw_response.content, "html.parser")
        search_results = soup.find_all("h4", class_="card-title")
        for i, result in enumerate(search_results[:3]):
            title = result.text.strip()
            function += f"Metabolomics Workbench Result {i+1}: {title}\n"
            link += f"Metabolomics Workbench Result {i+1}: [Insert link to the article here]\n"
            hit_count += 1

    # Search KEGG database for metabolite functions
    kegg_keywords = metabolite + ' ' + ' OR '.join(gut_keywords)
    kegg_search_url = f"https://www.kegg.jp/dbget-bin/www_bget?-f+m+compound+{kegg_keywords}"
    kegg_response = requests.get(kegg_search_url)
    if kegg_response.status_code == 200:
        soup = BeautifulSoup(kegg_response.content, "html.parser")
        kegg_result = soup.find("pre")
        if kegg_result:
            function += f"KEGG Result:\n{kegg_result.text.strip()}\n"
            link += f"KEGG Result: [Insert link to the article here]\n"
            hit_count += 1

    # Query PubChem API for metabolite functions
    pubchem_response = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{metabolite}/description/JSON")
    if pubchem_response.status_code == 200:
        pubchem_data = pubchem_response.json()
        if 'InformationList' in pubchem_data:
            pubchem_information = pubchem_data['InformationList']['Information']
            for info in pubchem_information:
                if 'Description' in info:
                    function += f"PubChem Result: {info['Description']}\n"
                if 'CID' in info:
                    link += f"PubChem Result: https://pubchem.ncbi.nlm.nih.gov/compound/{info['CID']}\n"
            hit_count += 1

    # Query PubMed for metabolite functions
    Entrez.email = "priyanshupandabcd@gmail.com"  # Replace with your email address
    pubmed_query = metabolite + ' ' + ' OR '.join(gut_keywords)
    handle = Entrez.esearch(db="pubmed", term=pubmed_query, retmax=3)
    pubmed_records = Entrez.read(handle)
    handle.close()
    pubmed_ids = pubmed_records["IdList"]

    for pubmed_id in pubmed_ids:
        pubmed_summary = None
        pubmed_link = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}"
        try:
            handle = Entrez.efetch(db="pubmed", id=pubmed_id, rettype="xml", retmode="text")
            record = Entrez.read(handle)
            pubmed_summary = record["PubmedArticle"][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0]
            handle.close()
        except:
            pass

        if pubmed_summary:
            if any(keyword.lower() in pubmed_summary.lower() for keyword in gut_keywords):
                function += f"PubMed Result: {pubmed_summary}\n"
                link += f"PubMed Result: {pubmed_link}\n"
                hit_count += 1

    # Search Google Scholar for metabolite functions related to the gut microbiome
    google_keywords = f"{metabolite} {' OR '.join(gut_keywords)}"
    google_search_url = f"https://scholar.google.com/scholar?q={google_keywords}"
    google_response = requests.get(google_search_url)
    if google_response.status_code == 200:
        soup = BeautifulSoup(google_response.content, "html.parser")
        search_results = soup.find_all("div", class_="gs_ri")
        for i, result in enumerate(search_results[:3]):
            title_element = result.find("h3", class_="gs_rt")
            abstract_element = result.find("div", class_="gs_rs")
            if title_element and abstract_element:
                title = title_element.text.strip()
                abstract = abstract_element.text.strip()
                if any(keyword.lower() in abstract.lower() for keyword in gut_keywords):
                    function += f"Google Scholar Result: {title}\n"
                    link += f"Google Scholar Result: [Insert link to the article here]\n"
                    hit_count += 1

    # Add metabolite, function, and hit count to respective lists
    functions.append(function)
    links.append(link)
    hits.append(hit_count)

    # If no gut keywords found, move to a separate sheet
    if hit_count == 0:
        no_keyword_df = no_keyword_df.append({'Metabolites': metabolite, 'Function': function}, ignore_index=True)

# Add new columns for function, links, and hits to the DataFrame
df['Function'] = functions
df['Links'] = links
df['Hits'] = hits

# Write the updated DataFrame back to the Excel file
with pd.ExcelWriter(excel_file) as writer:
    df.to_excel(writer, index=False, sheet_name='Main Sheet')
    no_keyword_df.to_excel(writer, index=False, sheet_name='No Keyword Mentions')
