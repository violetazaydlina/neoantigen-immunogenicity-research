import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import seaborn as sns

# Налаштування компактного академічного стилю
sns.set_theme(style="ticks")

# 1. Завантаження даних із файлів
df_netmhc = pd.read_csv("./data/processed/NECID_Results.csv")
df_neoad = pd.read_csv("./data/processed/foreigness score neoapred output.csv")

# Фільтруємо мутанти (беремо лише мутантні пептиди)
df_netmhc = df_netmhc[
    df_netmhc["peptide_orign"].astype(str).str.lower().str.contains("mut")
]

# Об'єднуємо таблиці за унікальним індексом рядка
merged_data = pd.merge(df_netmhc, df_neoad, on="necid_row_idx", how="inner")

# Виправляємо назву колонки response, якщо вона з суфіксом _x
resp_col = (
    "response_x"
    if "response_x" in merged_data.columns
    else ("response" if "response" in merged_data.columns else None)
)
if resp_col and resp_col != "response":
    merged_data = merged_data.rename(columns={resp_col: "response"})

# Шукаємо колонку з HLA автоматично
possible_hla_cols = [
    col
    for col in merged_data.columns
    if any(term in col.lower() for term in ["hla", "allele", "mhc"])
]
if not possible_hla_cols:
    raise ValueError("Не вдалося знайти колонку з HLA алелями в датасеті!")
hla_column = possible_hla_cols[0]

# 2. Збираємо дані для кожної алелі з урахуванням нових вимог
plot_rows = []

for hla_val, group in merged_data.groupby(hla_column):
    group_p = group[group["response"] == "P"]
    group_n = group[group["response"] == "N"]

    # Суворий фільтр: загалом >= 2 пептидів і обов'язково наявність і P, і N
    if len(group) < 2 or len(group_p) == 0 or len(group_n) == 0:
        continue

    # Розраховуємо співвідношення імуногенних до неімуногенних (P / N)
    ratio_p_n = len(group_p) / len(group_n)

    # Середнє значення Foreignness Score рахується по ВСІЙ групі алелі
    mean_ps = group["Foreignness_Score"].mean()

    plot_rows.append(
        {
            "HLA Allele": hla_val,
            "Mean_PS": mean_ps,
            "Ratio": ratio_p_n,
        }
    )

df_corr = pd.DataFrame(plot_rows)

if len(df_corr) > 2:
    # 3. Розрахунок коефіцієнта кореляції Спірмена між Mean_PS та співвідношенням P/N
    spearman_corr, p_value = stats.spearmanr(df_corr["Mean_PS"], df_corr["Ratio"])

    print("=" * 50)
    print(f"Коефіцієнт кореляції Спірмена (ρ): {spearman_corr:.4f}")
    print(f"p-value: {p_value:.4f}")
    print("=" * 50)

    # 4. Компактна наукова візуалізація
    plt.figure(figsize=(7.5, 4.5), dpi=300)

    # Точковий графік (Scatterplot)
    sns.scatterplot(
        data=df_corr,
        x="Mean_PS",
        y="Ratio",
        color="#2b8cbe",
        s=45,
        alpha=0.85,
        edgecolor="black",
        linewidth=0.5,
    )

    # Лінія тренду регресії
    sns.regplot(
        data=df_corr,
        x="Mean_PS",
        y="Ratio",
        scatter=False,
        color="#e34a33",
        ax=plt.gca(),
        line_kws={"linestyle": "--", "linewidth": 1.2},
    )

    # Акуратні підписи алелей біля кожної точки
    for _, row in df_corr.iterrows():
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

    # Налаштування осей, заголовків та сітки
    plt.title(
        f"Кореляція Спірмена: ρ = {spearman_corr:.2f} (p = {p_value:.3f})",
        fontsize=10,
        pad=10,
        weight="bold",
    )
    plt.xlabel("Середнє значення Foreignness Score ($ps$) в алелі", fontsize=9)
    plt.ylabel("Співвідношення (Імуногенні / Неімуногенні, $P/N$)", fontsize=9)

    plt.grid(True, linestyle=":", alpha=0.4, color="gray")
    sns.despine(top=True, right=True)
    
    plt.tight_layout()
    
    # За бажанням збереження у файл:
    # plt.savefig("hla_ratio_spearman_correlation.png", dpi=300, bbox_inches="tight")
    
    plt.show()
else:
    print("Занадто мало даних (менше 3 алелей) для розрахунку кореляції після фільтрації.")