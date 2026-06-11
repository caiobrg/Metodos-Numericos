import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Definição analítica da função objetivo para o plot das curvas de nível
def f_obj(x, y):
    return 2.0 * x * y + 2.0 * x - x**2 - 2.0 * y**2


def main():
    arquivos_necessarios = ["output1.dat", "output2.dat"]
    for arq in arquivos_necessarios:
        if not os.path.exists(arq):
            print(f"Erro: O arquivo '{arq}' não foi encontrado.")
            return

    print("Lendo os dados de otimização...")

    colunas_output = ["iter", "erro", "h", "x", "y", "dx", "dy"]
    df_am = pd.read_csv("output1.dat", sep=" ", header=None, names=colunas_output)
    df_gc = pd.read_csv("output2.dat", sep=" ", header=None, names=colunas_output)

    print("Gerando o gráfico...")
    fig, ax = plt.subplots(figsize=(10, 8))

    # Ponto ótimo conhecido
    x_opt, y_opt = 2.0, 1.0
    x0, y0 = df_am["x"].iloc[0], df_am["y"].iloc[0]

    # --- LIMITES EXATOS DO ESPAÇO ---
    # Pegamos o menor e maior valor de x e y que qualquer método alcançou
    min_x = min(df_am["x"].min(), df_gc["x"].min(), x_opt, x0)
    max_x = max(df_am["x"].max(), df_gc["x"].max(), x_opt, x0)
    min_y = min(df_am["y"].min(), df_gc["y"].min(), y_opt, y0)
    max_y = max(df_am["y"].max(), df_gc["y"].max(), y_opt, y0)

    # Margem de 20%
    margem_x = max((max_x - min_x) * 0.2, 1.0)
    margem_y = max((max_y - min_y) * 0.2, 1.0)

    lim_x = [min_x - margem_x, max_x + margem_x]
    lim_y = [min_y - margem_y, max_y + margem_y]

    # --- MALHA DA FUNÇÂO ---
    x_val = np.linspace(lim_x[0], lim_x[1], 400)
    y_val = np.linspace(lim_y[0], lim_y[1], 400)
    X, Y = np.meshgrid(x_val, y_val)
    Z = f_obj(X, Y)

    contorno = ax.contour(X, Y, Z, levels=40, alpha=0.8)
    ax.clabel(contorno, inline=True, fontsize=8, fmt="%.1f")

    # Plota as trajetórias
    ax.plot(
        df_am["x"],
        df_am["y"],
        marker="o",
        linestyle="-",
        markersize=4,
        linewidth=1.5,
        label="Aclive Máximo",
    )

    ax.plot(
        df_gc["x"],
        df_gc["y"],
        marker="s",
        linestyle="-",
        markersize=4,
        linewidth=1.5,
        label="Gradientes Conjugados (FR)",
    )

    # Marcação dos pontos principais
    ax.plot(x0, y0, marker="*", color="red", markersize=12, label="Ponto Inicial")
    ax.plot(
        x_opt,
        y_opt,
        marker="X",
        color="black",
        markersize=12,
        label="Ótimo Analítico (2, 1)",
    )

    # Ajusta os limites de vizualização para os pontos máximos e mnínimos
    ax.set_xlim(lim_x)
    ax.set_ylim(lim_y)

    # Configurações do gráfico
    ax.set_title(
        "Trajetórias de Otimização sobre as Curvas de Nível", fontsize=14, pad=15
    )
    ax.set_xlabel("Eixo X", fontsize=12)
    ax.set_ylabel("Eixo Y", fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.3)
    ax.legend(framealpha=0.9)

    plt.tight_layout()
    nome_saida = "trajetorias_otimizacao.png"
    plt.savefig(nome_saida, dpi=300)
    print(f"Gráfico salvo com sucesso como '{nome_saida}'.")
    plt.show()


if __name__ == "__main__":
    main()
