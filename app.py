import csv
import time
from io import StringIO
import gspread
import pandas as pd
import requests
from google.oauth2.service_account import Credentials


def predict_mhc_binding(sequence_text: str, allele: str) -> float:
    url = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
    seq_length = len(sequence_text.strip())

    payload = {
        "method": "netmhcpan-4.0",
        "sequence_text": sequence_text.strip(),
        "allele": f"HLA-{allele.strip()}",
        "length": str(seq_length),
    }

    response = requests.post(url, data=payload)
    response.raise_for_status()

    reader = csv.DictReader(StringIO(response.text.strip()), delimiter="\t")
    first_row = next(reader)

    return float(first_row["ic50"])


# 1. Авторизація
SCOPES = [
    "https://www.googleapis.com/auth/spreadsheets",
    "https://www.googleapis.com/auth/drive",
]
creds = Credentials.from_service_account_file("credentials.json", scopes=SCOPES)
client = gspread.authorize(creds)

# 2. Відкриття таблиці
spreadsheet = client.open_by_key("1DtZZqrSyvuir4e8MaTaJ6HM1gRCkka9xWt3HrPranyU")
sheet = spreadsheet.sheet1

records = sheet.get_all_records()

failed_items = []
recNum = 1  # Номер рядка в Google Sheets (заголовок у рядку 1, дані з рядка 2)

# ==========================================
# 1. ОСНОВНИЙ ЦИКЛ ОБРОБКИ
# ==========================================
for record in records:
    recNum += 1

    # Фільтрація записів
    if record.get("Imunogenicity") == "N":
        continue

    if record.get("Imunogenicity") == "P":
        if record.get("HLA II"):
            continue
        alele = record.get("HLA I")

    # ⚠️ Беремо Wild Type (якщо потрібен мутант — замініть на відповідне поле)
    minimalPeptide = record.get("Mutant peptide")

    if not minimalPeptide or not alele:
        continue

    print(f"Processing line #{recNum} (Peptide: {minimalPeptide}, Allele: {alele})")

    # Спроба №1 та №2
    try:
        result = predict_mhc_binding(sequence_text=minimalPeptide, allele=alele)
        sheet.update_acell(f"M{recNum}", result)
    except Exception as e:
        print(f"Mistake on line {recNum}: {e}. Retrying in 10s...")
        time.sleep(10)
        try:
            result = predict_mhc_binding(sequence_text=minimalPeptide, allele=alele)
            sheet.update_acell(f"M{recNum}", result)
        except Exception as e2:
            print(f"Failed again on line {recNum}. Saving to retry queue.")
            failed_items.append(
                {"row": recNum, "peptide": minimalPeptide, "allele": alele}
            )
if failed_items:
    print(f"\n--- Starting batch retry for {len(failed_items)} failed items ---")
    
    # Це офіційний API-ендпоінт IEDB, який містить ваш NetMHCpan-4.0 під капотом
    url = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
    
    fasta_lines = []
    alleles_set = set()
    row_mapping = {}  # Зв'язок: назва пептиду в пакеті -> номер рядка в Google Таблиці
    
    # Крок 1: Збираємо всі невдалі пептиди у форматі FASTA
    for idx, item in enumerate(failed_items, 1):
        p_name = f"peptide{idx}"
        
        # Додаємо пептид
        fasta_lines.append(f">{p_name}\n{item['peptide']}")
        
        # Форматуємо алель (гарантуємо наявність префікса HLA-)
        raw_allele = item['allele'].strip()
        formatted_allele = raw_allele if raw_allele.startswith("HLA-") else f"HLA-{raw_allele}"
        alleles_set.add(formatted_allele)
        
        # Запам'ятовуємо рядок таблиці
        row_mapping[p_name] = item['row']

    # Об'єднуємо все в єдиний текст для пакетного запиту
    sequence_text_batch = "\n".join(fasta_lines)
    alleles_batch = ",".join(alleles_set)

    # Параметри запиту (payload), де чітко вказано ваш метод netmhcpan-4.0
    payload = {
        "method": "netmhcpan-4.0",     # <--- Ваш метод прогнозування
        "sequence_text": sequence_text_batch,
        "allele": alleles_batch,
    }

    # Крок 2: Надсилаємо один спільний запит на всі помилкові елементи
    try:
        print("Sending batch request to IEDB server...")
        time.sleep(120)
        response = requests.post(url, data=payload)
        response.raise_for_status()

        # Крок 3: Читаємо відповідь і оновлюємо Google Таблицю
        reader = csv.DictReader(StringIO(response.text.strip()), delimiter="\t")
        
        for row in reader:
            seq_name = row.get("seq_name") or row.get("peptide_id")
            
            if seq_name in row_mapping:
                ic50_val = float(row["ic50"])
                target_row = row_mapping[seq_name]
                
                sheet.update_acell(f"M{target_row}", ic50_val)
                print(f"Successfully processed line #{target_row} via batch!")

    except Exception as e:
        print(f"Batch request failed: {e}")

print("\nDone!")
