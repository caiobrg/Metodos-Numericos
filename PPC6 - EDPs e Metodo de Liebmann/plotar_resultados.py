import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scienceplots
import os

plt.style.use("science")

# Parâmetros físicos da aleta conforme ConfigSimulacao
L = 0.5
H = 1.0
k = 200.0
h = 100.0
Tb = 425.0
Tinf = 273.0


def solucao_analitica(x):
    """Calcula a solução analítica 1D para a temperatura ao longo da aleta."""
    # Parâmetro da aleta 1D (perímetro/área = 2/H para o modelo 2D assumido)
    m = np.sqrt((2 * h) / (k * H))
    h_mk = h / (m * k)

    num_theta = np.cosh(m * (L - x)) + h_mk * np.sinh(m * (L - x))
    den_theta = np.cosh(m * L) + h_mk * np.sinh(m * L)

    theta = num_theta / den_theta
    T_ana = Tinf + (Tb - Tinf) * theta
    return T_ana


def plotar_relaxacao():
    """Lê os resultados de relaxação e plota iterações e tempo vs omega."""
    filename = "resultados/estudo_relaxacao.csv"
    if not os.path.exists(filename):
        print(f"Arquivo {filename} não encontrado.")
        return

    df = pl.read_csv(filename)

    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Eixo para iterações (vermelho)
    ax1.set_xlabel(r"Fator de Relaxação ($\omega$)")
    ax1.set_ylabel("Número de Iterações")
    sns.lineplot(
        data=df,
        x="omega",
        y="iteracoes",
        ax=ax1,
        marker="o",
        label="Iterações",
    )
    ax1.tick_params(axis="y")

    # Segundo Eixo para tempo (azul)
    ax2 = ax1.twinx()
    ax2.set_ylabel("Tempo (ms)")
    sns.lineplot(
        data=df,
        x="omega",
        y="tempo_ms",
        ax=ax2,
        marker="s",
        label="Tempo (ms)",
    )
    ax2.tick_params(axis="y")

    # Ajusta legendas para os dois eixos
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper center")
    if ax2.get_legend():
        ax2.get_legend().remove()

    plt.title(r"Efeito do Fator de Relaxação ($\omega$) na Convergência (Malha 81x81)")
    fig.tight_layout()
    plt.savefig("figures/grafico_relaxacao.png", dpi=300)
    print("Gráfico de relaxação salvo em 'figures/grafico_relaxacao.png'")
    # plt.show()
    plt.close("all")


def plotar_refinamento():
    """Lê os perfis numéricos e os compara com a solução analítica na linha central."""
    malhas = [11, 21, 41, 81]

    plt.figure(figsize=(10, 6))
    # Plot analítico contínuo
    x_ana = np.linspace(0, L, 200)
    T_ana = solucao_analitica(x_ana)
    plt.plot(x_ana, T_ana, "k--", label="Solução Analítica (1D)", linewidth=2)

    # Plot numérico para cada malha lendo os CSVs
    y_centro = H / 2.0

    for N in malhas:
        filename = f"resultados/temperatura_{N}x{N}.csv"
        if not os.path.exists(filename):
            print(f"Arquivo {filename} não encontrado.")
            continue

        df = pl.read_csv(filename)

        # Filtra os nós numéricos que estão na linha central y = 0.5 (com tolerância de precisão)
        df_centro = df.filter((pl.col("y") - y_centro).abs() < 1e-5).sort("x")

        x_num = df_centro["x"].to_numpy()
        T_num = df_centro["T"].to_numpy()

        plt.plot(x_num, T_num, marker=".", label=f"Numérico - Malha {N}x{N}")

    plt.title("Convergência de Malha vs Solução Analítica na Linha Central ($y=0.5m$)")
    plt.xlabel("Posição x (m)")
    plt.ylabel("Temperatura (K)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("figures/grafico_convergencia.png", dpi=300)
    print("Gráfico de convergência salvo em 'figures/grafico_convergencia.png'")
    # plt.show()
    plt.close("all")


def plotar_perfis_temperatura():
    """Plota os mapas de calor 2D para a distribuição de temperatura de cada malha."""
    malhas = [11, 21, 41, 81]

    for N in malhas:
        filename = f"resultados/temperatura_{N}x{N}.csv"
        if not os.path.exists(filename):
            continue

        df = pl.read_csv(filename)

        # Garante ordenação sequencial para o reshape matricial
        df = df.sort(["y", "x"])

        x = df["x"].to_numpy().reshape((N, N))
        y = df["y"].to_numpy().reshape((N, N))
        T = df["T"].to_numpy().reshape((N, N))

        plt.figure(figsize=(10, 5))
        # pcolormesh com shading='nearest' preserva a forma "quadriculada" real dos nós,
        # evidenciando totalmente a diferença visual de uma malha grossa para uma malha fina.
        contour = plt.pcolormesh(x, y, T, cmap="inferno", shading="nearest")
        plt.colorbar(contour, label="Temperatura (K)")

        plt.title(f"Perfil de Temperatura 2D - Malha {N}x{N}")
        plt.xlabel("Posição x (m)")
        plt.ylabel("Posição y (m)")
        plt.axis("equal")  # Mantém a proporção real da aleta (0.5 x 1.0)
        plt.tight_layout()

        out_filename = f"figures/perfil_temperatura_{N}x{N}.png"
        plt.savefig(out_filename, dpi=300)
        plt.close("all")
        print(f"Perfil de temperatura salvo em '{out_filename}'")


if __name__ == "__main__":
    if not os.path.exists("resultados"):
        print(
            "A pasta 'resultados' não existe. Certifique-se de executar o programa em C++ primeiro."
        )
    else:
        # Cria a pasta figures se não existir
        os.makedirs("figures", exist_ok=True)

        print("Gerando gráficos...")
        plotar_relaxacao()
        plotar_refinamento()
        plotar_perfis_temperatura()
