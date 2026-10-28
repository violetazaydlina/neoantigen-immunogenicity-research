import csv
import time
from io import StringIO

import pandas as pd
import requests

URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
TOOL_VERSION = "netmhcpan-4.1"
BATCH_SIZE = 10
PEP_COL = ["mut_peptide", "wt_peptide"]
immunogenicity = ["P", "N"]

# Tool constrains
unsupported_allelesA = ["A*24:01"]
min_peptide_length = 9
max_peptide_length = 11

print(f"Apply {TOOL_VERSION} to peptide={','.join(PEP_COL)}, responce={','.join(immunogenicity)}, peptide pairs length>={min_peptide_length}, peptide pairs length<={max_peptide_length}, HLA=alleleA, excluded alleleA: {','.join(unsupported_allelesA)} (unsuported)")

records = pd.read_csv("./data/raw/NECID_Query.csv")
#records[COL] = None

missing = records[
    (records["response"].isin(immunogenicity))
    & records["alleleA"].notna()
    & records["alleleB"].isna()
    & (~records["alleleA"].isin(unsupported_allelesA))
    & (
        (records["wt_peptide"].apply(len) <= max_peptide_length)
        & (records["wt_peptide"].apply(len) >= min_peptide_length)
    ) & (
        (records["mut_peptide"].apply(len) <= max_peptide_length)
        & (records["mut_peptide"].apply(len) >= min_peptide_length)
    )
]
missing = missing.drop_duplicates(PEP_COL+["alleleA"], keep="first")
print("MISSING", len(missing)*len(PEP_COL))
print(f"DETAILED")
print(missing["response"].value_counts()*len(PEP_COL))
print(missing["alleleA"].value_counts()*len(PEP_COL))
#print(missing["antigen_len"].value_counts()*len(PEP_COL))
for c in PEP_COL:
    print(missing[c].apply(len).value_counts())


# (row index, peptide, allele with HLA- prefix)
items = []
for i, row in missing.iterrows():
    allele = row["alleleA"].strip()
    allele = allele if allele.startswith("HLA-") else f"HLA-{allele}"
    immunogenicity = row["response"]
    for peptide_orign in PEP_COL:
        items.append((i, row[peptide_orign], allele, immunogenicity, peptide_orign))

results = []  # (row index, ic50)

for start in range(0, len(items), BATCH_SIZE):
    batch = items[start:start + BATCH_SIZE]
    print(f"Batch {start // BATCH_SIZE + 1}: {len(batch)} items")

    payload = {
        "method": TOOL_VERSION,
        "sequence_text": "\n".join(f">peptide{n}\n{pep}" for n, (_, pep, _, _, _) in enumerate(batch, 1)),
        "allele": ",".join(allele for _, _, allele, _, _ in batch),
        "length": ",".join(str(len(pep)) for _, pep, _, _, _ in batch),
    }

    #print(payload)

    try:
        time.sleep(15)  # be polite to the IEDB server
        response = requests.post(URL, data=payload)
        response.raise_for_status()
    except Exception as e:
        print(f"Batch failed: {e}")
        continue

    r = list(csv.DictReader(StringIO(response.text.strip()), delimiter="\t"))
    if len(r)<len(batch):
        print("RECIEVED")
        print(r)
        print("EXPECTED")
        print(payload)
        raise ValueError("query went wrong. see messege above.")
    else:
        # Lookup of predictions: (peptide, allele, length) -> ic50
        predictions = {
            (r["peptide"], r["allele"], int(r["length"])): float(r["ic50"])
            for r in csv.DictReader(StringIO(response.text.strip()), delimiter="\t")
        }

    for row_idx, pep, allele, immunogenicity, peptide_orign in batch:
        ic50 = predictions.get((pep, allele, len(pep)))
        if ic50 is not None:
            results.append([row_idx, immunogenicity, peptide_orign, allele, pep, ic50])
        else:
            print(f"No result for row {row_idx}: {pep} / {allele}")

print("WRITE RESULTS & SAVE", len(results))

results = pd.DtaFrame(results, columns = ["necid_row_idx", "response", "peptide_orign", "alleleA", "peptide_sequence", "ic50"])
results.to_csv("./data/processed/NECID_Results.csv", index=False)
print("Done!")