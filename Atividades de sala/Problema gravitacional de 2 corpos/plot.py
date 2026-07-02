import matplotlib.pyplot as plt
import numpy as np
import polars as pl

fig, axs = plt.subplots(ncols=2, figsize=(18, 6))

df = pl.read_csv("./results/resultados_euler.csv")

axs[0].plot(df["x1"], df["y1"], label="Corpo 1")
axs[0].plot(df["x2"], df["y2"], label="Corpo 2")
axs[0].set_title("Trajetórias")
axs[0].legend()

axs[1].plot(df["x1"] / df["x1"].max(), df["vx1"] / df["vx1"].max(), label="Corpo 1")
axs[1].plot(df["x2"] / df["x2"].max(), df["vx2"] / df["vx2"].max(), label="Corpo 2")
axs[1].set_title("Espaço de fase")
axs[1].legend()


fig.savefig("Resultados_Euler.png", bbox_inches="tight")
