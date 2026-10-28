import gspread
import matplotlib.pyplot as vialette
from google.oauth2.service_account import Credentials
import pandas as pd
import seaborn as sns

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

# Get all records as a list of dictionaries (using Row 1 as headers)
records = sheet.get_all_records()

data = []

recNum = 1

for record in records:
    recNum = recNum + 1
    print(f"processing line #{recNum}\n")
    if record['Imunogenicity'] == "N" or record['WT binding score'] == "":
        continue
    data.append(record['WT binding score'])

NEPdb = pd.DataFrame(records)
NEPdb["Binding score"] = pd.to_numeric(NEPdb["Binding score"], errors="coerce")
sns.boxplot(data=NEPdb, x="Imunogenicity", y="Binding score", hue="HLA")

vialette.boxplot(data)
vialette.title("Box Plot")
vialette.show()