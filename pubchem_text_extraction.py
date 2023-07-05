import requests
import xml.etree.ElementTree as ET
import pandas as pd

# Function to fetch the extracted text for a compound ID
def fetch_extracted_text(compound_id):
    # Base URL for PubChem compound data
    base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/XML/?response_type=display'

    # Format the URL with the compound ID
    url = base_url.format(compound_id)

    # Fetch the XML content
    response = requests.get(url)
    xml_content = response.content.decode('utf-8')

    # Parse the XML content
    root = ET.fromstring(xml_content)

    # Find the first five String elements
    string_elements = root.findall('.//{http://pubchem.ncbi.nlm.nih.gov/pug_view}String')[:15]
    extracted_text = [string_element.text.strip() for string_element in string_elements]
    return extracted_text 

# Read the Excel file
df = pd.read_excel(r'F:\NISER internship\Text-minning\pubchem_id_description\7179_metabolites\example.xlsx')

# Get the compound IDs from the second column
compound_ids = df.iloc[:, 1].tolist()

# Create a new column for extracted text
df['Extracted Text'] = ''

# Iterate over the compound IDs
for i, compound_id in enumerate(compound_ids, start=0):
    # Fetch the extracted text
    extracted_text = fetch_extracted_text(compound_id)
    
    # Store the extracted text in the corresponding row of the new column
    df.at[i, 'Extracted Text'] = extracted_text

# Save the modified DataFrame back to the Excel file
df.to_excel(r'F:\NISER internship\Text-minning\pubchem_id_description\7179_metabolites\example.xlsx', index=False)
