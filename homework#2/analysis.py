import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import re

df = pd.read_csv('homework#2/nextclade.csv', sep=';')
print(f"Loaded {len(df)} sequences")

# Variant distribution
print("\n" + "="*60)
print("Variant distrib")
print("="*60)
variant_counts = df['clade'].value_counts()
for clade, count in variant_counts.items():
    who = df[df['clade'] == clade].iloc[0]['clade_who']
    pct = (count / len(df) * 100)
    print(f"{clade} ({who}): {count} sequences ({pct:.1f}%)")

# Chronological? analysis
def extract_year(seqname):
    match = re.search(r'/(202[12])', seqname)
    return int(match.group(1)) if match else None

df['year'] = df['seqName'].apply(extract_year)

print("\n" + "="*60)
print("Chronological analysis")
print("="*60)
for year in [2021, 2022]:
    df_year = df[df['year'] == year]
    print(f"\n{year}: {len(df_year)} sequences")
    for clade in df_year['clade'].value_counts().index:
        count = len(df_year[df_year['clade'] == clade])
        pct = (count / len(df_year) * 100)
        who = df_year[df_year['clade'] == clade].iloc[0]['clade_who']
        print(f"  {clade} ({who}): {count} ({pct:.1f}%)")

# Spike mutations
print("\n" + "="*60)
print("spike m/tions")
print("="*60)
all_spike_muts = []
for aa_subs in df['aaSubstitutions']:
    if pd.notna(aa_subs):
        muts = [m for m in aa_subs.split(',') if m.startswith('S:')]
        all_spike_muts.extend(muts)

mut_counts = Counter(all_spike_muts)
spike_muts = pd.DataFrame({
    'Mutation': list(mut_counts.keys()),
    'Count': list(mut_counts.values()),
    'Percentage': [(c/len(df)*100) for c in mut_counts.values()]
}).sort_values('Count', ascending=False)

print("\n Highest Frequency Spike Mutations:")
print(spike_muts.head(15).to_string(index=False))

# Pie chart
fig, ax = plt.subplots(figsize=(12, 8))
colors = sns.color_palette("husl", len(variant_counts))

# labels = []
# for clade in variant_counts.index:
#     who = df[df['clade'] == clade].iloc[0]['clade_who']
#     labels.append(f"{clade} ({who})")

# wedges, texts, autotexts = ax.pie(
#     variant_counts, 
#     labels=labels,
#     autopct='%1.1f%%',
#     colors=colors,
#     startangle=90,
#     textprops={'fontsize': 11},
#     pctdistance=0.85
# )

# for text in texts:
#     text.set_fontsize(12)
#     text.set_fontweight('bold')

# for autotext in autotexts:
#     autotext.set_color('white')
#     autotext.set_fontsize(11)
#     autotext.set_fontweight('bold')

# ax.set_title('SARS-CoV-2 Variants in Kazakhstan\n(n=341 sequences)', 
#              fontsize=16, fontweight='bold', pad=20)

# plt.tight_layout()
# plt.savefig('homework#2/variant_distribution.png', dpi=300, bbox_inches='tight')
# print("\nSaved: variant_distribution.png")

# Bar chrrt
fig, ax = plt.subplots(figsize=(10, 6))
x_labels = [f"{c}\n({df[df['clade']==c].iloc[0]['clade_who']})" for c in variant_counts.index]
bars = ax.bar(range(len(variant_counts)), variant_counts.values, color=colors, edgecolor='black', linewidth=1.5)
ax.set_xticks(range(len(variant_counts)))
ax.set_xticklabels(x_labels, fontsize=11, fontweight='bold')
ax.set_ylabel('Number of Sequences', fontsize=13, fontweight='bold')
ax.set_title('SARS-CoV-2 Variants in Kazakhstan', fontsize=15, fontweight='bold', pad=15)
ax.grid(axis='y', alpha=0.3, linestyle='--')

for i, (bar, count) in enumerate(zip(bars, variant_counts.values)):
    height = bar.get_height()
    pct = (count / len(df) * 100)
    ax.text(bar.get_x() + bar.get_width()/2., height + 3,
            f'{int(count)}\n({pct:.1f}%)',
            ha='center', va='bottom', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig('homework#2/variant_distribution_bar.png', dpi=300, bbox_inches='tight')
print("Saved: variant_distribution_bar.png")

# Spike mutations
fig, ax = plt.subplots(figsize=(10, 8))
top_muts = spike_muts.head(15)
colors_vir = sns.color_palette("viridis", 15)
bars = ax.barh(range(len(top_muts)), top_muts['Count'], color=colors_vir, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top_muts)))
ax.set_yticklabels(top_muts['Mutation'], fontsize=10, fontweight='bold')
ax.invert_yaxis()
ax.set_xlabel('Number of Sequences', fontsize=12, fontweight='bold')
ax.set_title('Top 15 Spike Protein Mutations\n(n=341 sequences)', fontsize=14, fontweight='bold', pad=15)
ax.grid(axis='x', alpha=0.3, linestyle='--')

for i, (count, pct) in enumerate(zip(top_muts['Count'], top_muts['Percentage'])):
    ax.text(count + 5, i, f'{int(count)} ({pct:.1f}%)', 
            va='center', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.savefig('homework#2/spike_mutations.png', dpi=300, bbox_inches='tight')
print("Saved: spike_mutations.png")

# Chronological plot
fig, ax = plt.subplots(figsize=(10, 6))
years = ['2021\n(Aug-Dec)', '2022\n(Jan-May)']
df_2021 = df[df['year'] == 2021]
df_2022 = df[df['year'] == 2022]

all_clades = sorted(set(df['clade'].unique()))
bottom_2021 = 0
bottom_2022 = 0

for i, clade in enumerate(all_clades):
    count_2021 = len(df_2021[df_2021['clade'] == clade])
    count_2022 = len(df_2022[df_2022['clade'] == clade])
    
    color = colors[i] if i < len(colors) else 'gray'
    
    ax.bar(0, count_2021, bottom=bottom_2021, label=clade, color=color, edgecolor='black', linewidth=0.5)
    ax.bar(1, count_2022, bottom=bottom_2022, color=color, edgecolor='black', linewidth=0.5)
    
    bottom_2021 += count_2021
    bottom_2022 += count_2022

ax.set_xticks([0, 1])
ax.set_xticklabels(years, fontsize=12, fontweight='bold')
ax.set_ylabel('Number of Sequences', fontsize=12, fontweight='bold')
ax.set_title('Temporal Distribution of SARS-CoV-2 Variants', fontsize=14, fontweight='bold', pad=15)
ax.legend(title='Variant', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)

ax.text(0, len(df_2021) + 5, f'n={len(df_2021)}', ha='center', fontweight='bold', fontsize=11)
ax.text(1, len(df_2022) + 5, f'n={len(df_2022)}', ha='center', fontweight='bold', fontsize=11)

plt.tight_layout()
plt.savefig('homework#2/temporal_distribution.png', dpi=300, bbox_inches='tight')
print("Saved: temporal_distribution.png")

spike_muts.to_csv('homework#2/spike_mutations.csv', index=False)
print("Saved: spike_mutations.csv")

print("\n" + "="*60)
print("done~!")
print("="*60)
print("\ncreated files")
print("  - barchgrt.png")
print("  - chrnlgy.png")
print("  - spike_mutations.png")
print("  - tableofspikemutations.csv")