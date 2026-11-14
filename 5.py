import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Налаштування компактного академічного стилю
sns.set_theme(style="ticks")

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

# 2. Збираємо співвідношення імуногенних (P) до неімуногенних (N) по кожній алелі
plot_rows = []

for hla_val, group in merged_data.groupby(hla_column):
    group_p = group[group["response"] == "P"]
    group_n = group[group["response"] == "N"]

    if len(group) < 2 or len(group_p) == 0 or len(group_n) == 0:
        continue

    count_p = len(group_p)
    count_n = len(group_n)

    # Розраховуємо співвідношення (Імуногенні / Неімуногенні)
    if count_n > 0:
        ratio = count_p / count_n
    else:
        ratio = count_p  # Якщо неімуногенних немає, беремо кількість P

    # Середнє ps беремо загалом по всій алелі (або за бажанням можна окремо для P/N)
    mean_ps = group["Foreignness_Score"].mean()

    if count_p > 0 or count_n > 0:
        plot_rows.append(
            {
                "HLA Allele": hla_val,
                "Mean_PS": mean_ps,
                "Ratio": ratio,
                "Count_P": count_p,
                "Count_N": count_n,
            }
        )

df_plot = pd.DataFrame(plot_rows)

if df_plot.empty:
    print("Увага: Немає даних для побудови.")
else:
    # 3. Компактне полотно для зручного перегляду та збереження
    plt.figure(figsize=(7.5, 4.5), dpi=300)

    # Зменшені акуратні точки
    sns.scatterplot(
        data=df_plot,
        x="Mean_PS",
        y="Ratio",
        color="#2b8cbe",
        s=45,
        alpha=0.85,
        edgecolor="black",
        linewidth=0.5,
    )

    # Лінія тренду (регресії)
    if len(df_plot) >= 2:
        sns.regplot(
            data=df_plot,
            x="Mean_PS",
            y="Ratio",
            scatter=False,
            color="#e34a33",
            ax=plt.gca(),
            line_kws={"linestyle": "--", "linewidth": 1.2},
        )

    # Компактні підписи алелей меншим шрифтом (fontsize=6)
    for _, row in df_plot.iterrows():
        plt.annotate(
            row["HLA Allele"],
            (row["Mean_PS"], row["Ratio"]),
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

    plt.title(
        "Залежність співвідношення (Імуногенні / Неімуногенні) від Foreignness Score ($ps$)",
        fontsize=10,
        pad=10,
        weight="bold",
    )
    plt.xlabel("Середнє значення Foreignness Score ($ps$) в алелі", fontsize=9)
    plt.ylabel("Співвідношення (Імуногенні / Неімуногенні)", fontsize=9)

    plt.grid(True, linestyle=":", alpha=0.4, color="gray")
    sns.despine(top=True, right=True)

    plt.tight_layout()
    plt.show()
