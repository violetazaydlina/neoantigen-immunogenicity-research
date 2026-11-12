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


# 4. Підрахунок кількості для NECID (позитивні та негативні)
nec_pos = (merged_data["response"] == "P").sum()  # Імуногенні (NECID)
nec_neg = (merged_data["response"] == "N").sum()  # Неімуногенні (NECID)

# Створюємо датасет для побудови
plot_data = pd.DataFrame(
    {
        "Category": ["Імуногенні (NECID)", "Неімуногенні (NECID)"],
        "Count": [nec_pos, nec_neg],
    }
)

# 5. Побудова графіка
plt.figure(figsize=(8, 6))
sns.barplot(data=plot_data, x="Category", y="Count", palette="Set2")

plt.title("Розподіл імуногенності за базою NECID")
plt.xlabel("Категорія")
plt.ylabel("Кількість пептидів")
plt.tight_layout()

# Виведення та збереження
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

# Створюємо датасет для нового деталізованого графіка
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

# 5 (продовження). Побудова згрупованого графіка
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

print(
    f"Імуногенні з високим Foreignness (>0.5): {pos_high}"
) 
print(
    f"Імуногенні з низьким Foreignness (<=0.5): {pos_low}"
)  
print(
    f"Неімуногенні з високим Foreignness (>0.5): {neg_high}"
) 
print(
    f"Неімуногенні з низьким Foreignness (<=0.5): {neg_low}"
)  

# Виведення та збереження другого графіка
plt.show()
