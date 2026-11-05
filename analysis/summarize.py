"""Summarize the supplied CSVs and draw descriptive figures; no MD required."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/summary'
FIG = ROOT / 'figures'
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'svg.fonttype': 'none'})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / f'{name}.png', dpi=180, bbox_inches='tight')
    plt.close(fig)


def main():
    metrics = pd.read_csv(ROOT / 'data/boltz/metrics.csv')
    labels = [f'Binder{i}' for i in range(1, 6)] + ['Brazzein']
    keys = ['confidence_score', 'iptm', 'complex_plddt', 'complex_pde']
    stats = metrics.groupby('label')[keys].agg(['count', 'mean', 'std'])
    stats.columns = ['_'.join(c) for c in stats.columns]
    stats.reindex(labels).to_csv(OUT / 'boltz_summary.csv')
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.6))
    for ax, key in zip(axes, keys):
        for x, name in enumerate(labels):
            values = metrics.loc[metrics.label == name, key].to_numpy()
            ax.scatter(x + np.linspace(-0.11, 0.11, len(values)), values, s=20, alpha=.7)
            ax.plot([x-.2, x+.2], [values.mean()]*2, color='black', lw=2)
        ax.set_xticks(range(6), labels, rotation=45, ha='right')
        ax.set_title(key)
        ax.set_ylabel('Predicted distance error (Å)' if key.endswith('pde') else 'Confidence (0–1)')
    fig.suptitle('Boltz-1: five model predictions per complex; black line = mean', y=1.04)
    save(fig, 'boltz_summary')

    energy_rows = []
    energy_values = []
    systems = ['binder1', 'binder2', 'binder3', 'advantame']
    for system in systems:
        df = pd.read_csv(ROOT / f'data/energies/{system}_mmgbsa.csv')
        values = df['r_psp_MMGBSA_dG_Bind']
        assert len(values) == 90 and values.notna().all()
        energy_values.append(values.to_numpy())
        energy_rows.append(dict(system=system, n_results=len(values),
                                mean_kcal_mol=values.mean(), sample_sd_kcal_mol=values.std(ddof=1),
                                median_kcal_mol=values.median(), intended_window_ns='20–200',
                                note='Snapshot timestamps not present in the result CSV'))
    pd.DataFrame(energy_rows).to_csv(OUT / 'mmgbsa_summary.csv', index=False)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.boxplot(energy_values, tick_labels=['Binder1', 'Binder2', 'Binder3', 'Advantame'], showfliers=False)
    rng = np.random.default_rng(42)
    for x, values in enumerate(energy_values, 1):
        ax.scatter(x+rng.uniform(-.12, .12, len(values)), values, alpha=.35, s=12)
        ax.scatter(x, values.mean(), marker='D', color='black', s=25, zorder=5)
    ax.set_ylabel('Prime MM/GBSA score (kcal/mol)')
    ax.set_title('90 results per system; diamond = mean\nBoxes: median and quartiles; whiskers: 1.5 × IQR')
    save(fig, 'mmgbsa_summary')

    fig, axes = plt.subplots(3, 3, figsize=(12, 8))
    for i in range(1, 4):
        for j, key in enumerate(['rg', 'rmsd', 'rmsf']):
            df = pd.read_csv(ROOT / f'data/md/binder{i}/{key}.csv')
            ax = axes[i-1, j]
            y = 'RMSF CA (A)' if key == 'rmsf' else df.columns[1]
            ax.plot(df.iloc[:, 0], df[y], lw=.8)
            ax.set_title(f'Binder{i}: {key.upper()}')
            ax.set_xlabel('Residue number' if key == 'rmsf' else 'Time (ns)')
            ax.set_ylabel('Å')
    save(fig, 'binder_md_summary')

    ranges = []
    summaries = []
    for path in sorted((ROOT / 'data/md').glob('*/*.csv')):
        df = pd.read_csv(path)
        time_col = next((c for c in df.columns if c.lower().startswith('time')), None)
        if time_col is None:
            continue
        time = pd.to_numeric(df[time_col])
        ranges.append(dict(file=str(path.relative_to(ROOT)), rows=len(df),
                           start_ns=time.min(), end_ns=time.max()))
        for start in [0, 20]:
            window = df[(time >= start-.01) & (time <= 200.01)]
            for col in window.select_dtypes(include='number'):
                if col == time_col or col == 'frame':
                    continue
                values = window[col].dropna()
                summaries.append(dict(file=str(path.relative_to(ROOT)), metric=col,
                                      requested_window_ns=f'{start}–200',
                                      available_start_ns=window[time_col].min(),
                                      available_end_ns=window[time_col].max(),
                                      n=len(values), mean=values.mean(),
                                      sample_sd=values.std(ddof=1)))
    pd.DataFrame(ranges).to_csv(OUT / 'time_ranges.csv', index=False)
    pd.DataFrame(summaries).to_csv(OUT / 'md_descriptive_summary.csv', index=False)
    print(f'Wrote descriptive summaries to {OUT} and three overview figures to {FIG}')


if __name__ == '__main__':
    main()
