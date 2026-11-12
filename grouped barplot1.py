import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df_netmhc = pd.read_csv("./data/processed/NECID_Results.csv")
# Додаємо шлях ./data/processed/ або перевірте назву файлу
df_neoad = pd.read_csv("./data/processed/foreigness score neoapred output.csv")

merged_data = pd.merge(df_netmhc, df_neoad, on="necid_row_idx", how="inner")

resp_col = (
    "response_x"
    if "response_x" in merged_data.columns
    else ("response" if "response" in merged_data.columns else None)
)

counts = (
    merged_data[["peptide_orign", resp_col]]
    .value_counts()
    .to_frame()
    .reset_index()
)
counts.columns = ["peptide_orign", "response", "count"]

sns.barplot(data=counts, x="peptide_orign", y="count", hue="response")

plt.show()