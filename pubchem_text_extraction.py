import requests
import xml.etree.ElementTree as ET
import pandas as pd

# URL template for PubChem XML API
url_template = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/XML/?response_type=display'

# Read the Excel file
excel_file = "F:/NISER internship/Text-minning/pubchem_id_description/metabolite_name_1.xlsx"
df = pd.read_excel(excel_file)

# Iterate over the rows of the DataFrame
for index, row in df.iterrows():
    # Get the PubChem ID from the 'pubchem_id' column
    pubchem_id = row['PubChem ID']

    # Create the URL for the specific PubChem ID
    url = url_template.format(pubchem_id)

    # Fetch the XML content                                                     
    response = requests.get(url)
    xml_content = response.content.decode('utf-8')

    # Parse the XML content
    root = ET.fromstring(xml_content)

    # Find the first five String elements
    string_elements = root.findall('.//{http://pubchem.ncbi.nlm.nih.gov/pug_view}String')[:5]

    # Extract the text under each String element
    string_texts = [element.text.strip() for element in string_elements]

    # Join the extracted texts with a delimiter
    string_texts_combined = ' '.join(string_texts)

    # Update the 'function' column in the DataFrame
    df.at[index, 'function'] = string_texts_combined

# Save the updated DataFrame back to the Excel file
df.to_excel(excel_file, index=False)
