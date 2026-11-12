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
# Перевірте, чи назва ключа збігається у ваших файлах (тут використано 'necid_row_idx')
merged_df = pd.merge(
    df_netmhc[df_netmhc["peptide_orign"] == "mut_peptide"],  # беремо мутантні пептиди
    df_neoad,
    on="necid_row_idx",
    how="inner",
    suffixes=("_netmhc", "_neoad"),
)

print(f"Успішно об'єднано даних. Загальна кількість пар: {len(merged_df)}")

# 4. Очищення та підготовка даних для аналізу
# Фільтруємо невалідні значення та перетворюємо на числовий формат
merged_df["ic50"] = pd.to_numeric(merged_df["ic50"], errors="coerce")
merged_df["foreignness_score"] = pd.to_numeric(
    merged_df["foreignness_score"], errors="coerce"
)
merged_df = merged_df.dropna(subset=["ic50", "foreignness_score", "response"])

# Для зручності візуалізації IC50 часто беруть у логарифмічній шкалі або обмежують викиди
# (сильні зв'язувачі мають низький IC50, тому логарифм або інверсія часто покращують сприйняття)
merged_df["log_ic50"] = merged_df["ic50"].apply(lambda x: max(x, 0.1))

# 5. Створення комплексного графіку (2 панелі / підграфіки)
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# --- Панель А: Scatter Plot (Foreignness Score vs IC50 з розподілом на імунний статус) ---
sns.scatterplot(
    data=merged_df,
    x="foreignness_score",
    y="ic50",
    hue="response",
    palette={"P": "#2ca02c", "N": "#d62728"},  # Зелений для P, Червоний для N
    alpha=0.7,
    s=60,
    ax=axes[0],
)

# Додаємо лінії порогових значень (наприклад, IC50 = 50 нМ для сильних зв'язувачів)
axes[0].axhline(
    y=50,
    color="gray",
    linestyle="--",
    alpha=0.7,
    label="Поріг сильної афінності (50 нМ)",
)
axes[0].set_yscale("log")  # Логарифмічна шкала для IC50
axes[0].set_title(
    "А) Зіставлення індексу чужорідності та афінності ($\mathbf{IC_{50}}$)",
    fontsize=13,
    pad=10,
)
axes[0].set_xlabel("Просторовий індекс чужорідності (Foreignness Score)", fontsize=11)
axes[0].set_ylabel("Афінність зв'язування $\mathbf{IC_{50}}$ (нМ, лог. шкала)", fontsize=11)
axes[0].legend(title="Імунний статус", labels=["Позитивний (P)", "Негативний (N)"])

# --- Панель Б: Boxplot / Violinplot (Розподіл Foreignness Score залежно від імунного статусу) ---
sns.boxplot(
    data=merged_df,
    x="response",
    y="foreignness_score",
    palette=["#2ca02c", "#d62728"],
    ax=axes[1],
    width=0.4,
    boxprops=dict(alpha=0.8),
)
sns.stripplot(
    data=merged_df,
    x="response",
    y="foreignness_score",
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

# Збереження зображення у високій якості для наукової роботи
output_path = "./data/processed/Neoantigen_Comprehensive_Analysis.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"Графік успішно зберечено у файл: {output_path}")

# Відображення графіка
plt.show()
input()