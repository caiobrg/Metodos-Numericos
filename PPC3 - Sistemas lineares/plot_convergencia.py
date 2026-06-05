import argparse

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Plota o estudo de convergencia log-log"
    )
    parser.add_argument(
        "--arq_esp", type=str, required=True, help="CSV de convergencia espacial"
    )
    parser.add_argument(
        "--arq_temp", type=str, required=True, help="CSV de convergencia temporal"
    )
    parser.add_argument("--L", type=float, required=True, help="Espessura L")
    args = parser.parse_args()

    # ==========================================
    # 1. LER DADOS DE CONVERGÊNCIA ESPACIAL
    # ==========================================
    df_esp = pd.read_csv(args.arq_esp)
    N_nos = df_esp["N"].values
    R2_esp = df_esp["R2"].values

    dx = args.L / (N_nos - 1)
    erro_esp = np.abs(1.0 - R2_esp)

    # ==========================================
    # 2. LER DADOS DE CONVERGÊNCIA TEMPORAL
    # ==========================================
    df_temp = pd.read_csv(args.arq_temp)
    dt = df_temp["dt"].values
    R2_temp = df_temp["R2"].values

    erro_temp = np.abs(1.0 - R2_temp)

    # ==========================================
    # PLOTAGEM ESPACIAL
    # ==========================================
    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.loglog(dx, erro_esp, "bo-", linewidth=2, markersize=6, label="Erro Simulado")

    # Calculando a ordem do método espacial (usando apenas os 3 primeiros pontos)
    if len(dx) >= 3:
        log_dx = np.log(dx[0:3])
        log_err_esp = np.log(erro_esp[0:3])
        m_esp, c_esp = np.polyfit(log_dx, log_err_esp, 1)
        plt.loglog(
            dx,
            np.exp(m_esp * np.log(dx) + c_esp),
            "k--",
            alpha=0.7,
            label=rf"Ordem calculada: $O(\Delta x^{{{m_esp:.2f}}})$",
        )

    plt.title(r"Convergência Espacial ($\Delta t$ fixo)", fontweight="bold")
    plt.xlabel(r"Tamanho da Malha ($\Delta x$)")
    plt.ylabel("Erro ($1 - R^2$)")
    plt.grid(True, which="both", ls=":", alpha=0.6)
    plt.legend()
    plt.gca().invert_xaxis()

    # ==========================================
    # PLOTAGEM TEMPORAL
    # ==========================================
    plt.subplot(1, 2, 2)
    plt.loglog(dt, erro_temp, "ro-", linewidth=2, markersize=6, label="Erro Simulado")

    # Calculando a ordem do método temporal (isolando a zona estável dt <= 0.125)
    indices_estaveis = np.where(dt <= 0.125)[0]
    if len(indices_estaveis) >= 2:
        log_dt = np.log(dt[indices_estaveis])
        log_err_temp = np.log(erro_temp[indices_estaveis])
        m_temp, c_temp = np.polyfit(log_dt, log_err_temp, 1)
        plt.loglog(
            dt,
            np.exp(m_temp * np.log(dt) + c_temp),
            "k--",
            alpha=0.7,
            label=rf"Ordem calculada: $O(\Delta t^{{{m_temp:.2f}}})$",
        )

    plt.title("Convergência Temporal ($N$ fixo)", fontweight="bold")
    plt.xlabel(r"Passo de Tempo ($\Delta t$)")
    plt.grid(True, which="both", ls=":", alpha=0.6)
    plt.legend()
    plt.gca().invert_xaxis()

    plt.tight_layout()
    plt.savefig("estudo_convergencia.png", dpi=300)
    print("Gráficos de convergência salvos em 'estudo_convergencia.png'.")
    # plt.show()


if __name__ == "__main__":
    main()
