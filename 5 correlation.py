import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import seaborn as sns

# Налаштування академічного стилю
sns.set_theme(style="ticks")

# 1. Завантаження даних із файлів
df_netmhc = pd.read_csv("./data/processed/NECID_Results.csv")
df_neoad = pd.read_csv("./data/processed/foreigness score neoapred output.csv")

# ВИПРАВЛЕНО: Правильна фільтрація мутантів за допомогою str.contains
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

plot_rows = []

for hla_val, group in merged_data.groupby(hla_column):
    # 1. Формуємо підгрупи класів
    group_p = group[group["response"] == "P"]
    group_n = group[group["response"] == "N"]

    # 2. ЖОРСТКИЙ ФІЛЬТР: якщо умови не виконуються, одразу пропускаємо цю алель
    if len(group) < 2 or len(group_p) == 0 or len(group_n) == 0:
        continue

    # 3. Якщо пройшло фільтр — додаємо дані для подальшого аналізу та графіка
    plot_rows.append(
        {
            "HLA Allele": hla_val,
            "Mean_PS": group_p["Foreignness_Score"].mean(),
            "Count": len(group_p),
        }
    )

df_corr = pd.DataFrame(plot_rows)

if len(df_corr) > 2:
    # 3. Розрахунок кореляції за Спірменом
    spearman_corr, p_value = stats.spearmanr(df_corr["Mean_PS"], df_corr["Count"])

    print("=" * 50)
    print(f"Коефіцієнт кореляції Спірмена (ρ): {spearman_corr:.4f}")
    print(f"p-value: {p_value:.4f}")
    print("=" * 50)

    # 4. Компактна візуалізація
    plt.figure(figsize=(7, 4.5), dpi=300)

    sns.scatterplot(
        data=df_corr,
        x="Mean_PS",
        y="Count",
        color="#d95f02",
        s=40,
        alpha=0.85,
        edgecolor="black",
        linewidth=0.5,
    )

    # Підписи алелей
    for _, row in df_corr.iterrows():
        plt.annotate(
            row["HLA Allele"],
            (row["Mean_PS"], row["Count"]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=6,
            alpha=0.8,
            bbox=dict(
                boxstyle="round,pad=0.1",
                fc="white",
                ec="lightgray",
                alpha=0.6,
                linewidth=0.4,
            ),
        )

    # Виводимо значення ρ на графік
    plt.title(
        f"Кореляція Спірмена: ρ = {spearman_corr:.2f} (p = {p_value:.3f})",
        fontsize=10,
        pad=10,
        weight="bold",
    )
    plt.xlabel("Середнє значення Foreignness Score ($ps$)", fontsize=9)
    plt.ylabel("Кількість імуногенних пептидів ($n$)", fontsize=9)

    plt.grid(True, linestyle=":", alpha=0.4, color="gray")
    sns.despine(top=True, right=True)
    plt.tight_layout()
    plt.show()
else:
    print("Занадто мало даних для розрахунку кореляції після фільтрації.")