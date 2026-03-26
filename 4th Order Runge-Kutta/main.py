import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re

plt.rcParams["figure.dpi"] = 400
plt.rcParams["figure.figsize"] = (10, 6)
plt.rcParams["figure.autolayout"] = True
sns.set_theme(style="whitegrid", palette="viridis", context="paper")


def analytic_solution(t, St, Re):
    if Re == 0.0:
        return 1.0 - np.exp(-t / St)

    v_star = (-1.0 + np.sqrt(1.0 + 1.5 * Re)) / (0.75 * Re)
    Q = (3.0 / (8.0 * St)) * Re
    P = -(3.0 / (4.0 * St)) * Re * v_star - (1.0 / St)

    exp_term = np.exp(-P * t)
    inverse_term = (Q / P) + ((-1.0 / v_star) - (Q / P)) * exp_term
    v_z = v_star + (1.0 / inverse_term)

    return v_z


def plot_stokes_comparison(directory="results/St", figureDir="results/St"):
    dir_path = Path(directory)
    fig_dir_path = Path(figureDir)

    # Searches for all the .csv files in the directory
    csv_files = list(dir_path.glob("*.csv"))

    if not csv_files:
        print(
            f"Warning: no .csv file found in specified directory '{directory}'."
        )
        return

    plt.figure()

    for file in csv_files:
        # Regex to extract the St value from file name
        match = re.search(r"St_([0-9.]+)\.csv", file.name)
        if not match:
            print(f"Ignoring {file.name}: Name pattern not recognised.")
            continue

        st_val = float(match.group(1))

        # Importing the data
        df = pd.read_csv(file)
        t_num = df["Time(t*)"].to_numpy()
        v_num = df["Velocity(v*)"].to_numpy()

        v_ana = analytic_solution(t=t_num, St=st_val, Re=0.0)

        # Maximum absolute error (measure for the difference between analytic and numerical solutions)
        err_max = np.max(np.abs(v_num - v_ana))

        # Figure plot
        _ = sns.lineplot(
            x=t_num,
            y=v_ana,
            label=f"Analítica (St={st_val})",
            linewidth=2,
            zorder=1,
        )
        _ = sns.lineplot(
            x=t_num,
            y=v_num,
            linestyle="dashed",
            label=f"Numérica (St={st_val} | Erro Max={err_max:.2e})",
        )

    _ = plt.title(
        "Relaxação da Velocidade (Regime de Stokes)", fontsize=14, pad=15
    )
    _ = plt.xlabel("Tempo Adimensional (t*)", fontsize=12)
    _ = plt.ylabel("Velocidade Adimensional (v*)", fontsize=12)
    _ = plt.legend(frameon=True, facecolor="white", framealpha=0.9)
    _ = plt.tight_layout()

    # Salva a imagem gerada dentro da própria subpasta
    plt.savefig(fig_dir_path / "comparação_St.png")
    plt.close("all")


def plot_dt_analysis_completa(
    directory="results/Dt", figureDir="results/Dt", St=1.0, Re=0.5
):
    dir_path = Path(directory)
    fig_dir_path = Path(figureDir)
    csv_files = list(dir_path.glob("*.csv"))

    if not csv_files:
        print(f"Warning: no .csv file found in '{directory}'.")
        return

    # Regex function to extract the DT value from file name
    def extract_dt(f: Path) -> float:
        m = re.search(r"Dt_([0-9.]+)\.csv", f.name, re.IGNORECASE)
        return float(m.group(1)) if m else -1.0

    # Ordering the file names (to look better  on the legend and to make code logic easier)
    csv_files_sorted = sorted(
        [f for f in csv_files if extract_dt(f) > 0],
        key=extract_dt,
        reverse=True,
    )

    err_list = []
    dts_list = []

    # 2 figures, one for comparison between the curves and other for convergence analysis
    fig_conv, ax_conv = plt.subplots()
    fig_comp, ax_comp = plt.subplots()

    t_num = np.array([])
    v_ana = np.array([])

    # Importing data (for both the figures)
    for file in csv_files_sorted[::-1]:
        dt_val = extract_dt(file)
        df = pd.read_csv(file)

        t_num = df["Time(t*)"].to_numpy()
        v_num = df["Velocity(v*)"].to_numpy()

        # Error measure, same as in Stokes comparison
        v_ana = analytic_solution(t_num, St, Re)
        erro_max = np.max(np.abs(v_num - v_ana))

        dts_list.append(dt_val)
        err_list.append(erro_max)

        # Plots on the comparison graph
        _ = sns.lineplot(
            x=t_num,
            y=v_num,
            linestyle="--",
            alpha=0.8,
            label=f"Numérica (dt={dt_val})",
            zorder=2,
            ax=ax_comp,
        )

    # Final adjustments on graph 1
    # Last "t_num" and "v_ana" are used as they have the highest amount of points
    _ = sns.lineplot(
        x=t_num,
        y=v_ana,
        label=f"Analítica (St={St}, Re={Re})",
        linewidth=2,
        zorder=1,
        color="black",
        ax=ax_comp,
    )
    _ = ax_comp.set_title(
        "Efeito do Passo de Tempo (dt) na Curva de Relaxação",
        fontsize=14,
        pad=15,
    )
    ax_comp.set_xlabel("Tempo Adimensional (t*)", fontsize=12)
    ax_comp.set_ylabel("Velocidade Adimensional (v*)", fontsize=12)
    ax_comp.legend(frameon=True, facecolor="white", framealpha=0.9)
    fig_comp.tight_layout()
    fig_comp.savefig(fig_dir_path / "comparacao_curvas_dt.png", dpi=300)
    ## End of graph 1 (comparison) ##

    # plotting on graph 2
    dts_list = np.array(dts_list)
    err_list = np.array(err_list)

    _ = sns.scatterplot(
        x=dts_list,
        y=err_list,
        color="red",
        s=100,
        label="Erro Numérico (RK4)",
        zorder=3,
        ax=ax_conv,
    )
    _ = sns.lineplot(
        x=dts_list, y=err_list, color="red", alpha=0.5, zorder=2, ax=ax_conv
    )

    anchor_cte = err_list[-1] / (dts_list[-1] ** 4)
    dt4_ref = anchor_cte * (dts_list**4)  ## dt^4 reference line

    _ = sns.lineplot(
        x=dts_list,
        y=dt4_ref,
        color="gray",
        linestyle="dashed",
        label=r"Referência Teórica $\mathcal{O}(dt^4)$",
        zorder=1,
        ax=ax_conv,
    )

    ax_conv.set_xscale("log")
    ax_conv.set_yscale("log")
    ax_conv.set_title(
        "Convergência do RK4 com Efeito Inercial", fontsize=14, pad=15
    )
    ax_conv.set_xlabel("Passo de Tempo Adimensional (dt)", fontsize=12)
    ax_conv.set_ylabel(r"Erro Global Máximo (Norma $L_\infty$)", fontsize=12)
    ax_conv.grid(True, which="both", ls="--", alpha=0.4)
    ax_conv.legend(frameon=True, facecolor="white", framealpha=0.9)
    fig_conv.tight_layout()
    fig_conv.savefig(fig_dir_path / "analise_convergencia_dt.png")
    ## End of graph 2 (convergence) ##

    # Closes all figures as to not interfere in future figures
    plt.close("all")


def plot_reynolds_comparison(
    directory="results/Re", figureDir="results/Re", st_fixo=1.0
):
    dir_path = Path(directory)
    fig_dir_path = Path(figureDir)
    csv_files = list(dir_path.glob("*.csv"))

    if not csv_files:
        print(f"Warning: no .csv file found in '{directory}'.")
        return

    plt.figure()

    for file in csv_files:
        # Extrai o Re do nome do arquivo (ex: resultados_Re_0.5.csv)
        match = re.search(r"Re_([0-9.]+)\.csv", file.name, re.IGNORECASE)
        if not match:
            continue

        re_val = float(match.group(1))

        df = pd.read_csv(file)
        t_num = df["Time(t*)"].to_numpy()
        v_num = df["Velocity(v*)"].to_numpy()

        # Solução Analítica Exata para o Re da iteração atual
        v_ana = analytic_solution(t_num, st_fixo, re_val)

        err_max = np.max(np.abs(v_num - v_ana))

        _ = sns.lineplot(
            x=t_num,
            y=v_ana,
            label=f"Analítica (Re={re_val})",
            linewidth=2,
        )
        _ = sns.lineplot(
            x=t_num,
            y=v_num,
            linestyle="dashed",
            alpha=1,
            label=f"Numérica (Re={re_val} | Erro Max={err_max:.2e})",
        )

    plt.title(
        "Relaxação da Velocidade (Efeito Inercial - Arrasto de Oseen)",
        fontsize=14,
        pad=15,
    )
    plt.xlabel("Tempo Adimensional (t*)", fontsize=12)
    plt.ylabel("Velocidade Adimensional (v*)", fontsize=12)
    plt.legend(frameon=True, facecolor="white", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(fig_dir_path / "comparacao_Re.png")
    plt.close("all")


plot_stokes_comparison("results/St", "documento/figures")
plot_dt_analysis_completa("results/Dt", "documento/figures")
plot_reynolds_comparison("results/Re", "documento/figures")
