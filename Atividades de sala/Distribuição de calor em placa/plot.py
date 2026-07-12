import os

import matplotlib.pyplot as plt
import pandas as pd

for file in os.listdir("results/"):
    df = pd.read_csv("results/" + file)

    # Transforma a lista 1D em uma Matriz 2D (Pivot)
    malha_matriz = df.pivot(index="y", columns="x", values="T")

    plt.figure(figsize=(10, 6))

    # 1. Desenha o fundo colorido (Substituindo o heatmap do Seaborn)
    # origin="lower" garante que o Y=0 comece na base do gráfico automaticamente
    # extent define os limites reais do gráfico nos eixos [X_min, X_max, Y_min, Y_max]
    mapa_fundo = plt.imshow(
        malha_matriz.values,
        cmap="coolwarm",
        origin="lower",
        extent=[  # pyright: ignore[reportArgumentType]
            malha_matriz.columns.min(),
            malha_matriz.columns.max(),
            malha_matriz.index.min(),
            malha_matriz.index.max(),
        ],
        aspect="auto",
    )

    # Adiciona a barra lateral de cores (Legenda de temperatura)
    plt.colorbar(mapa_fundo, label="Temperatura [K]")

    # 2. Adiciona as Isotermas (Agora no mesmo sistema de coordenadas!)
    cs = plt.contour(
        malha_matriz.columns,
        malha_matriz.index,
        malha_matriz.values,
        levels=13,
        colors="black",
        linewidths=1.2,
    )

    # 3. Rótulos das linhas de contorno
    plt.clabel(cs, inline=True, fmt="%1.1f K", fontsize=9)

    plt.title(f"Campo de Temperatura e Isotermas - {file}")
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")

plt.show()
