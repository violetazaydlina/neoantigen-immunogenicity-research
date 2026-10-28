import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

data = pd.read_csv("./data/processed/NECID_Results.csv")

counts = data[["peptide_orign", "response"]].value_counts().to_frame().reset_index()
sns.barplot(data=counts, x="peptide_orign", y="count", hue="response");
plt.savefig("./plots/necid_barplot_oregin_respose.png", format="png")
plt.close()

sns.boxplot(data=data, x="peptide_orign", y="ic50", hue="response");
plt.savefig("./plots/necid_boxplot_oregin_respose_ic50.png", format="png")
plt.close()

def get_binder_category(ic50):
    ic50 = float(ic50)
    if ic50<=50:
        return "strong-binder"
    elif ic50 >50 and ic50<=500:
        return "weak-binder"
    elif ic50>500:
        return "non-binder"
    
data["binder_category"] = data["ic50"].apply(get_binder_category)
COL="mut_peptide"
counts = data[data["peptide_orign"]==COL][["binder_category", "response"]].value_counts().to_frame().reset_index()
sns.barplot(data=counts, x="binder_category", y="count", hue="response", order=["strong-binder","weak-binder","non-binder"]);
plt.title(COL)
plt.savefig(f"./plots/necid_barplot_binder_category_respose_{COL}.png", format="png")
plt.close()

COL="mut_peptide"
counts = data[data["peptide_orign"]==COL][["alleleA", "response"]].value_counts().to_frame().reset_index()
plt.figure(figsize=(12,6))
sns.barplot(data=counts, x="alleleA", y="count", hue="response", order=sorted(counts['alleleA'].unique()));
plt.title(COL)
plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
plt.tight_layout()
plt.savefig(f"./plots/necid_barplot_allele_respose_{COL}.png", format="png")
plt.close()

plt.figure(figsize=(12,6))
sns.boxplot(data=data[data["peptide_orign"]==COL], x="alleleA", y="ic50", hue="response", order=sorted(data[data["peptide_orign"]==COL]['alleleA'].unique()));
plt.title(COL)
plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
plt.tight_layout()
plt.savefig(f"./plots/necid_boxplot_allele_respose_ic50_{COL}.png", format="png")
plt.close()

counts = data[data["peptide_orign"]==COL][["binder_category","alleleA", "response"]].value_counts().to_frame().reset_index()
# what allele has the most strong-binders in response=N than response=P
# percent of strong-binders per allele and response
counts["allele-binder-percent"] = round((
    counts["count"] / counts.groupby(["alleleA", "response"])["count"].transform("sum")
)*100, 2)
COL_2 = "strong-binder"
plt.figure(figsize=(12,6))
sns.barplot(data=counts[counts["binder_category"]==COL_2], x="alleleA", y="allele-binder-percent", hue="response", order=sorted(counts['alleleA'].unique()));
plt.title(f"{COL}, {COL_2}")
plt.ylabel("percent")
plt.xticks(rotation=45, ha='right', rotation_mode='anchor')
plt.tight_layout()
plt.savefig(f"./plots/necid_barplot_allele_respose_{COL}_{COL_2}.png", format="png")
plt.close()