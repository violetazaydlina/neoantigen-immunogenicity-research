import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# 1. Налаштування стилю візуалізації (для наукових робіт)
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["font.size"] = 11

# 2. Завантаження обох датасетів
print("Завантаження даних...")
df_netmhc = pd.read_csv("./data/processed/NECID_Results.csv")
df_neoad = pd.read_csv("./data/processed/foreigness score neoapred output.csv")

# 3. Об'єднання таблиць (злиття за унікальним індексом зразка / ID)
merged_df = pd.merge(
    df_netmhc[
        df_netmhc["peptide_orign"].astype(str).str.lower().str.contains("mut")
    ],
    df_neoad,
    on="necid_row_idx",
    how="inner",
    suffixes=("_netmhc", "_neoad"),
)

print(f"Успішно об'єднано даних. Загальна кількість пар: {len(merged_df)}")

# 4. Очищення та підготовка даних для аналізу
merged_df["ic50"] = pd.to_numeric(merged_df["ic50"], errors="coerce")
merged_df["foreignness_score"] = pd.to_numeric(
    merged_df["foreignness_score"], errors="coerce"
)
merged_df = merged_df.dropna(subset=["ic50", "foreignness_score", "response"])

# 5. Створення комплексного графіку (2 панелі / підграфіки)
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Загальна єдина палітра для всього скрипта: Зелений для P (позитивні), Червоний для N (негативні)
palette_dict = {"P": "#2ca02c", "N": "#d62728"}

# --- Панель А: Scatter Plot (Foreignness Score vs IC50) ---
sns.scatterplot(
    data=merged_df,
    x="foreignness_score",
    y="ic50",
    hue="response",
    hue_order=["P", "N"],  # Суворо фіксуємо порядок для відповідності легенді
    palette=palette_dict,
    alpha=0.7,
    s=60,
    ax=axes[0],
)

axes[0].axhline(
    y=50,
    color="gray",
    linestyle="--",
    alpha=0.7,
    label="Поріг сильної афінності (50 нМ)",
)
axes[0].set_yscale("log")
axes[0].set_title(
    "А) Зіставлення індексу чужорідності та афінності ($\mathbf{IC_{50}}$)",
    fontsize=13,
    pad=10,
)
axes[0].set_xlabel("Просторовий індекс чужорідності (Foreignness Score)", fontsize=11)
axes[0].set_ylabel("Афінність зв'язування $\mathbf{IC_{50}}$ (нМ, лог. шкала)", fontsize=11)

# Виправляємо порядок елементів у легенді панелі А
handles, labels = axes[0].get_legend_handles_labels()
# Зберігаємо лінію порогу (якщо вона є) і коректно підписуємо P та N
axes[0].legend(
    handles=handles,
    labels=[
        "Поріг сильної афінності (50 нМ)",
        "Позитивні неоантигени (P)",
        "Негативні пептиди (N)",
    ],
    title="Статус та поріг",
)

# --- Панель Б: Boxplot / Stripplot (Розподіл Foreignness Score за статусом) ---
sns.boxplot(
    data=merged_df,
    x="response",
    y="foreignness_score",
    order=["P", "N"],  # Спочатку P (зелений), потім N (червоний)
    palette=palette_dict,
    ax=axes[1],
    width=0.4,
    boxprops=dict(alpha=0.8),
)
sns.stripplot(
    data=merged_df,
    x="response",
    y="foreignness_score",
    order=["P", "N"],
    color="black",
    alpha=0.3,
    jitter=0.2,
    size=4,
    ax=axes[1],
)

axes[1].set_title(
    "Б) Розподіл індексу чужорідності за імунним статусом", fontsize=13, pad=10
)
axes[1].set_xlabel("Експериментальний імунний статус", fontsize=11)
axes[1].set_xticklabels(
    ["Позитивні неоантигени (P)", "Негативні пептиди (N)"]
)
axes[1].set_ylabel("Просторовий індекс чужорідності (Foreignness Score)", fontsize=11)

# Загальне компонування
plt.tight_layout()

# Збереження зображення
output_path = "./data/processed/Neoantigen_Comprehensive_Analysis.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"Графік успішно зберечено у файл: {output_path}")

plt.show()