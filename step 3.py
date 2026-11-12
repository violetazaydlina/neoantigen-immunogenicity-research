import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# 1. Завантаження даних
df_netmhc = pd.read_csv("./data/processed/NECID_Results.csv")
df_neoad = pd.read_csv("./data/processed/foreigness score neoapred output.csv")

df_netmhc = df_netmhc[
    df_netmhc["peptide_orign"].astype(str).str.lower().str.contains("mut")
]

print(f"Залишилося рядків після фільтрації мутантів: {len(df_netmhc)}")

# 2. Об'єднання за індексом
merged_data = pd.merge(df_netmhc, df_neoad, on="necid_row_idx", how="inner")
print(f"Кількість рядків після об'єднання з прогнозами: {len(merged_data)}")

# 3. Визначення правильної колонки response
resp_col = (
    "response_x"
    if "response_x" in merged_data.columns
    else ("response" if "response" in merged_data.columns else None)
)
if resp_col and resp_col != "response":
    merged_data = merged_data.rename(columns={resp_col: "response"})


# 4. Підрахунок кількості для NECID загалом (позитивні та негативні)
nec_pos = (merged_data["response"] == "P").sum()
nec_neg = (merged_data["response"] == "N").sum()

plot_data = pd.DataFrame(
    {
        "Category": ["Імуногенні (NECID)", "Неімуногенні (NECID)"],
        "Count": [nec_pos, nec_neg],
    }
)

# Побудова першого загального графіка (розподіл NECID)
plt.figure(figsize=(8, 6))
sns.barplot(data=plot_data, x="Category", y="Count", palette="Set2")
plt.title("Розподіл імуногенності за базою NECID (Загальний)")
plt.xlabel("Категорія")
plt.ylabel("Кількість пептидів")
plt.tight_layout()
plt.show()


pos_high = (
    (merged_data["response"] == "P") & (merged_data["Foreignness_Score"] > 0.5)
).sum()
pos_low = (
    (merged_data["response"] == "P") & (merged_data["Foreignness_Score"] <= 0.5)
).sum()
neg_high = (
    (merged_data["response"] == "N") & (merged_data["Foreignness_Score"] > 0.5)
).sum()
neg_low = (
    (merged_data["response"] == "N") & (merged_data["Foreignness_Score"] <= 0.5)
).sum()

plot_data_detailed = pd.DataFrame(
    {
        "NECID Status": [
            "Імуногенні (P)",
            "Імуногенні (P)",
            "Неімуногенні (N)",
            "Неімуногенні (N)",
        ],
        "Foreignness Group": [
            "> 0.5 (High)",
            "<= 0.5 (Low)",
            "> 0.5 (High)",
            "<= 0.5 (Low)",
        ],
        "Count": [pos_high, pos_low, neg_high, neg_low],
    }
)

# Другий загальний деталізований графік
plt.figure(figsize=(10, 6))
sns.barplot(
    data=plot_data_detailed,
    x="NECID Status",
    y="Count",
    hue="Foreignness Group",
    palette="Set2",
)
plt.title(
    "Детальний розподіл пептидів за NECID та порогом Foreignness Score (0.5)"
)
plt.xlabel("Експериментальний статус (NECID)")
plt.ylabel("Кількість пептидів")
plt.legend(title="Foreignness Score")
plt.tight_layout()
plt.show()

print(f"Імуногенні з високим Foreignness (>0.5) [Загалом]: {pos_high}")
print(f"Імуногенні з низьким Foreignness (<=0.5) [Загалом]: {pos_low}")
print(f"Неімуногенні з високим Foreignness (>0.5) [Загалом]: {neg_high}")
print(f"Неімуногенні з низьким Foreignness (<=0.5) [Загалом]: {neg_low}")


# =====================================================================
# ЗБІР І ПІДСУМКОВИЙ ГРАФІК У РОЗРІЗІ ВСІХ HLA АЛЕЛЬ
# =====================================================================
print("\n--- Збір статистики та побудова підсумкового графіка за HLA ---")

possible_hla_cols = [
    col
    for col in merged_data.columns
    if any(term in col.lower() for term in ["hla", "allele", "mhc"])
]

if possible_hla_cols:
    hla_column = possible_hla_cols[0]
    print(f"Знайдено колонку для HLA: '{hla_column}'")

    summary_rows = []

    for hla_val, group in merged_data.groupby(hla_column):
        h_pos_high = (
            (group["response"] == "P") & (group["Foreignness_Score"] > 0.5)
        ).sum()
        h_pos_low = (
            (group["response"] == "P") & (group["Foreignness_Score"] <= 0.5)
        ).sum()
        h_neg_high = (
            (group["response"] == "N") & (group["Foreignness_Score"] > 0.5)
        ).sum()
        h_neg_low = (
            (group["response"] == "N") & (group["Foreignness_Score"] <= 0.5)
        ).sum()

        # Друкуємо числа в термінал для кожного HLA (без зайвих графіків)
        print(
            f"[HLA: {hla_val}] TP(Pos High): {h_pos_high}, FN(Pos Low): {h_pos_low} | FP(Neg High): {h_neg_high}, TN(Neg Low): {h_neg_low}"
        )

        # Додаємо результати до загального списку для підсумкового графіка
        summary_rows.append(
            {
                "HLA": hla_val,
                "Category": "Імуногенні (>0.5)",
                "Count": h_pos_high,
            }
        )
        summary_rows.append(
            {"HLA": hla_val, "Category": "Імуногенні (<=0.5)", "Count": h_pos_low}
        )
        summary_rows.append(
            {
                "HLA": hla_val,
                "Category": "Неімуногенні (>0.5)",
                "Count": h_neg_high,
            }
        )
        summary_rows.append(
            {
                "HLA": hla_val,
                "Category": "Неімуногенні (<=0.5)",
                "Count": h_neg_low,
            }
        )

    # Робимо один спільний підсумковий датафрейм
    df_summary = pd.DataFrame(summary_rows)

    # Побудова ЄДИНОГО підсумкового графіка для всіх HLA у кінці
    plt.figure(figsize=(14, 7))
    sns.barplot(
        data=df_summary,
        x="HLA",
        y="Count",
        hue="Category",
        palette="viridis",
    )
    plt.title("Порівняльний розподіл пептидів за типами HLA та категоріями")
    plt.xlabel("Тип HLA алелі")
    plt.ylabel("Кількість пептидів")
    plt.xticks(rotation=45, ha="right")
    plt.legend(title="Категорії", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()

else:
    print(
        "Увага: Не вдалося автоматично знайти колонку з назвою HLA/allele/mhc."
    )
