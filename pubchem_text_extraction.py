import requests
import xml.etree.ElementTree as ET
import pandas as pd

# API key
api_key = "62815109555cb6e1c776695ca04e03cae608"

# URL template for PubChem XML API
url_template = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{{}}/XML/?response_type=display&api_key={api_key}'

# Read the Excel file
excel_file = "F:/NISER internship/Text-minning/pubchem_id_description/New folder/metabolite_name_1.xlsx"
df = pd.read_excel(excel_file)

# Create empty list to store the extracted text
extracted_text = []

# Iterate over the rows of the DataFrame
for index, row in df.iterrows():
    # Get the PubChem ID from the 'PubChem ID' column
    pubchem_id = row['PubChem ID']

    # Create the URL for the specific PubChem ID with the API key
    url = url_template.format(pubchem_id)

    # Fetch the XML content
    response = requests.get(url)
    xml_content = response.content.decode('ISO-8859-1')

    # Parse the XML content
    root = ET.fromstring(xml_content)

    # Find the first five String elements
    string_elements = root.findall('.//{http://pubchem.ncbi.nlm.nih.gov/pug_view}String')[:30]

    # Extract the text from the String elements
    strings = [string_element.text.strip() for string_element in string_elements]

    # Append the extracted text to the list
    extracted_text.append(strings)

# Add the extracted text as a new column in the DataFrame
df['Extracted_Text'] = extracted_text

# Save the updated DataFrame back to the Excel file
df.to_excel(excel_file, index=False)
