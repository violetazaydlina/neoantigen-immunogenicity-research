import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# 1. Завантаження даних із файлів
df_netmhc = pd.read_csv("./data/processed/NECID_Results.csv")
df_neoad = pd.read_csv("./data/processed/foreigness score neoapred output.csv")

# Фільтруємо мутанти
df_netmhc = df_netmhc[
    df_netmhc["peptide_orign"].astype(str).str.lower().str.contains("mut")
]
merged_data = pd.merge(df_netmhc, df_neoad, on="necid_row_idx", how="inner")

# Виправляємо назву колонки response, якщо вона з суфіксом _x
resp_col = (
    "response_x"
    if "response_x" in merged_data.columns
    else ("response" if "response" in merged_data.columns else None)
)
if resp_col and resp_col != "response":
    merged_data = merged_data.rename(columns={resp_col: "response"})

# Шукаємо колонку з HLA
possible_hla_cols = [
    col
    for col in merged_data.columns
    if any(term in col.lower() for term in ["hla", "allele", "mhc"])
]
hla_column = possible_hla_cols[0]

# 2. Збираємо метрики для кожної алелі з урахуванням фільтру (>= 10 пептидів)
heatmap_rows = []

for hla_val, group in merged_data.groupby(hla_column):
    total_peptides = len(group)


    group_p = group[group["response"] == "P"]
    group_n = group[group["response"] == "N"]
    if len(group) < 2 or len(group_p) == 0 or len(group_n) == 0:
           continue


    
    # Матриця помилок
    tp = ((group["response"] == "P") & (group["Foreignness_Score"] > 0.5)).sum()
    fn = ((group["response"] == "P") & (group["Foreignness_Score"] <= 0.5)).sum()
    fp = ((group["response"] == "N") & (group["Foreignness_Score"] > 0.5)).sum()
    tn = ((group["response"] == "N") & (group["Foreignness_Score"] <= 0.5)).sum()

    # Метрики якості
    accuracy = (tp + tn) / total_peptides
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = (
        (2 * precision * recall) / (precision + recall)
        if (precision + recall) > 0
        else 0
    )

    immunogenic_count = tp + fn
    non_immunogenic_count = fp + tn

    ratio_im_nonim = (
        immunogenic_count / non_immunogenic_count
        if non_immunogenic_count > 0
        else immunogenic_count
    )

    heatmap_rows.append(
        {
            "HLA Allele": hla_val,
            "Accuracy": accuracy,
            "Precision": precision,
            "Recall": recall,
            "F1 Score": f1,
            "Im/Non-Im Ratio": ratio_im_nonim,
            "Total Neoantigens": total_peptides,
        }
    )

# 3. Перетворюємо у DataFrame
df_heatmap = pd.DataFrame(heatmap_rows)

if df_heatmap.empty:
    print(
        "Увага: Немає жодної алелі, яка б містила 10 або більше пептидів."
    )
else:
    df_heatmap.set_index("HLA Allele", inplace=True)

    # Нормалізуємо кожну колонку самостійно від 0 до 1 для коректного фарбування кольорів
    df_normalized = (df_heatmap - df_heatmap.min()) / (
        df_heatmap.max() - df_heatmap.min()
    )
    # Якщо якийсь стовпчик має однакові значення (ділення на 0), заповнюємо його нулями
    df_normalized = df_normalized.fillna(0)

    # 4. Побудова теплової карти (Heatmap)
    plt.figure(figsize=(10, max(6, len(df_heatmap) * 0.5)))
    sns.heatmap(
        df_normalized,
        annot=df_heatmap,  # Виводимо оригінальні реальні числа всередині комірок!
        fmt=".2f",  # Формат виведення чисел
        cmap="YlGnBu",  # Колірна палитра
        linewidths=0.5,
        cbar_kws={"label": "Нормалізований рівень (0 до 1)"},
    )

    plt.title(
        "Теплова карта метрик у розрізі HLA алелей (n >= 10)",
        fontsize=13,
        pad=15,
    )
    plt.xlabel("Метрики та показники", fontsize=11)
    plt.ylabel("Група HLA алель", fontsize=11)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()

    plt.show()
