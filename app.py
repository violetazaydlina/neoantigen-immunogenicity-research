import gspread
import requests
import csv
import time
from io import StringIO

from google.oauth2.service_account import Credentials

def predict_mhc_binding(sequence_text: str, allele: str) -> float:
    """Sends a request to the IEDB MHC-I binding prediction API.

    Args:
        sequence_text: The amino acid sequence (e.g., "ARFTGIKTA").
        allele: The target allele (e.g., "HLA-A*01:01").

    Returns:
        The raw text response from the IEDB API.
    """
    url = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"

    # Automatically calculate sequence length
    seq_length = len(sequence_text.strip())

    payload = {
        "method": "netmhcpan-4.0",
        "sequence_text": sequence_text.strip(),
        "allele": f"HLA-{allele.strip()}",
        "length": str(seq_length),
    }

    response = requests.post(url, data=payload)

    # Raise an exception if the request was unsuccessful (e.g., 4xx, 5xx)
    response.raise_for_status()

    # Parse response
    reader = csv.DictReader(StringIO(response.text.strip()), delimiter="\t")
    first_row = next(reader)

    return float(first_row["ic50"])


# 1. Define the scope
SCOPES = [
    "https://www.googleapis.com/auth/spreadsheets",
    "https://www.googleapis.com/auth/drive",
]

# 2. Authenticate using the service account JSON
creds = Credentials.from_service_account_file("credentials.json", scopes=SCOPES)
client = gspread.authorize(creds)

# 3. Open the spreadsheet by name or by URL/Key
# Option A: By Name
# spreadsheet = client.open("My Spreadsheet Name")

# Option B: By Key (found in the sheet's URL)
spreadsheet = client.open_by_key("1DtZZqrSyvuir4e8MaTaJ6HM1gRCkka9xWt3HrPranyU")

# Select the active worksheet (e.g., first sheet)
sheet = spreadsheet.sheet1


# --- READING DATA ---

# Get all values as a 2D list of strings
# all_values = sheet.get_all_values()
# print("All Data (Lists):", all_values)

# Get all records as a list of dictionaries (using Row 1 as headers)
records = sheet.get_all_records()
# print("All Data (Dicts):", records)

recNum = 1

for record in records:
    recNum = recNum + 1
    print(f"processing line #{recNum}\n")
    if record['Imunogenicity'] == "N":
        continue
    #if recNum < 574: continue
    # if recNum == 25:
    #     break
    # print (record['Minimal peptide'])
    minimalPeptide = record['WT peptide']
    
    if record['Imunogenicity'] == "P":
        if (record['HLA II']):
            continue
        alele = record['HLA I']
        try:
            result = predict_mhc_binding(sequence_text=minimalPeptide, allele=alele)
        except Exception as e:
            print(f"error: {e}")
            time.sleep(5)
            try:
                result = predict_mhc_binding(sequence_text=minimalPeptide, allele=alele)
            except: 
                print(f"error: {e}")
                continue            
            
    
        # Update a single cell (e.g., cell B2)
        sheet.update_acell(f"N{recNum}", result)
    else:
        continue

    

