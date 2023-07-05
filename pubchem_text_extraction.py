import requests
import xml.etree.ElementTree as ET
import pandas as pd

# Function to fetch the extracted text for a list of compound IDs
def fetch_extracted_texts(compound_ids):
    extracted_texts = []
    # Iterate over the compound IDs
    for compound_id in compound_ids:
        # Base URL for PubChem compound data
        base_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{}/XML/?response_type=display'

        # Format the URL with the compound ID
        url = base_url.format(compound_id)

        # Fetch the XML content
        response = requests.get(url)
        xml_content = response.content.decode('latin-1')  # Decode using 'latin-1'

        # Parse the XML content
        root = ET.fromstring(xml_content)

        # Find the first five String elements
        string_elements = root.findall('.//{http://pubchem.ncbi.nlm.nih.gov/pug_view}String')[:100]
        extracted_text = [string_element.text.strip() for string_element in string_elements]
        extracted_texts.append(extracted_text)
    
    return extracted_texts

# Read the Excel file
df = pd.read_excel(r'F:\NISER internship\Text-minning\pubchem_id_description\7179_metabolites\New folder\Metabolites_with_pubchem_ID(5068).xlsx')

# Get the compound IDs from the second column
compound_ids = df.iloc[:, 1].tolist()

# Create a new column for extracted text
df['Extracted Text'] = ''

# Set the batch size
batch_size = 100

# Process compound IDs in batches
for i in range(0, len(compound_ids), batch_size):
    # Get the current batch of compound IDs
    batch_compound_ids = compound_ids[i:i+batch_size]

    # Fetch the extracted texts for the current batch
    extracted_texts = fetch_extracted_texts(batch_compound_ids)

    # Update the extracted texts in the corresponding rows of the DataFrame
    for j, extracted_text in enumerate(extracted_texts):
        df.at[i+j, 'Extracted Text'] = extracted_text

    # Save the modified DataFrame to a new Excel file after processing each batch
    df.to_excel(r'F:\NISER internship\Text-minning\pubchem_id_description\7179_metabolites\New folder\Metabolites_with_pubchem_ID_processed.xlsx', index=False)

# Save the final modified DataFrame back to the original Excel file
df.to_excel(r'F:\NISER internship\Text-minning\pubchem_id_description\7179_metabolites\New folder\Metabolites_with_pubchem_ID_processed.xlsx', index=False)
