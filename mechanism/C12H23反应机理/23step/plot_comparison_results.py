#!/usr/bin/env python3
"""
Plot comparison results for the three mechanisms:
- Three-body Pure (Original)
- Three-body Fixed  
- Cleaned (Remove FORD)
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

# Results from the latest comparison
phi_values = [0.8, 0.9, 1.0, 1.1, 1.2]

# Flame speed results (m/s)
results = {
    'Three-body Pure (Original)': [None, None, None, None, None],  # All failed
    'Three-body Fixed': [0.0044, 0.0057, 0.0071, 0.0085, 0.0098],
    'Cleaned (Remove FORD)': [0.0044, 0.0057, 0.0071, 0.0085, 0.0098]
}

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Flame speed comparison
ax1.plot(phi_values, results['Three-body Fixed'], 'bo-', linewidth=2, markersize=8, label='Three-body Fixed')
ax1.plot(phi_values, results['Cleaned (Remove FORD)'], 'rs-', linewidth=2, markersize=8, label='Cleaned (Remove FORD)')

# Mark failed cases for Three-body Pure
ax1.scatter(phi_values, [0.01]*len(phi_values), c='red', marker='x', s=100, label='Three-body Pure (Original) - Failed')

ax1.set_xlabel('Equivalence Ratio (φ)', fontsize=12)
ax1.set_ylabel('Laminar Flame Speed (m/s)', fontsize=12)
ax1.set_title('Flame Speed Comparison: C12H23 Mechanisms', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)
ax1.set_xlim(0.75, 1.25)

# Plot 2: Relative difference (showing identical results for working mechanisms)
working_mechanisms = ['Three-body Fixed', 'Cleaned (Remove FORD)']
for mech in working_mechanisms:
    ax2.plot(phi_values, results[mech], 'o-', linewidth=2, markersize=8, label=mech)

ax2.set_xlabel('Equivalence Ratio (φ)', fontsize=12)
ax2.set_ylabel('Laminar Flame Speed (m/s)', fontsize=12)
ax2.set_title('Working Mechanisms Detail View', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=10)
ax2.set_xlim(0.75, 1.25)

plt.tight_layout()
plt.savefig('flame_speed_comparison_updated.png', dpi=300, bbox_inches='tight')
plt.savefig('flame_speed_comparison_updated.pdf', bbox_inches='tight')

# Create a summary table figure
fig2, ax = plt.subplots(figsize=(12, 4))
ax.axis('tight')
ax.axis('off')

# Create table data
table_data = []
for phi in phi_values:
    row = [f"{phi:.1f}"]
    for mech in ['Three-body Pure (Original)', 'Three-body Fixed', 'Cleaned (Remove FORD)']:
        val = results[mech][phi_values.index(phi)]
        if val is None:
            row.append("Failed")
        else:
            row.append(f"{val:.4f}")
    table_data.append(row)

# Create table
col_labels = ['φ', 'Three-body Pure (Original)', 'Three-body Fixed', 'Cleaned (Remove FORD)']
table = ax.table(cellText=table_data, colLabels=col_labels, cellLoc='center', loc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.5)

# Color coding
colors = ['lightcoral' if cell == 'Failed' else 'lightgreen' for row in table_data for cell in row[1:]]
for i, (key, cell) in enumerate(table.get_celld().items()):
    if key[0] > 0 and key[1] > 0:  # Skip header row and first column
        if table_data[key[0]-1][key[1]] == 'Failed':
            cell.set_facecolor('lightcoral')
        else:
            cell.set_facecolor('lightgreen')

ax.set_title('Flame Speed Results Summary (m/s)', fontsize=14, fontweight='bold', pad=20)
plt.savefig('results_table.png', dpi=300, bbox_inches='tight')

# Create convergence analysis plot
fig3, ax = plt.subplots(figsize=(10, 6))

# Plot showing the convergence issue
mechanisms_test = ['Three-body Pure (Original)', 'Three-body Fixed', 'Cleaned (Remove FORD)']
convergence_status = [0, 1, 1]  # 0=failed, 1=success

bars = ax.bar(mechanisms_test, convergence_status, color=['red', 'green', 'green'], alpha=0.7)
ax.set_ylabel('Convergence Status (1=Success, 0=Failed)', fontsize=12)
ax.set_title('Mechanism Convergence Analysis', fontsize=14, fontweight='bold')
ax.set_ylim(0, 1.2)

# Add labels on bars
for bar, status in zip(bars, convergence_status):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
            'Success' if status == 1 else 'Failed',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig('convergence_analysis.png', dpi=300, bbox_inches='tight')

print("All plots have been generated successfully!")
print("Generated files:")
print("- flame_speed_comparison_updated.png/pdf")
print("- results_table.png")
print("- convergence_analysis.png")