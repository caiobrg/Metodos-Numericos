import os

import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import scienceplots  # noqa: F401

plt.style.use("science")

fig, axs = plt.subplots(ncols=3, figsize=(30, 6))

for file in os.listdir("results"):
    df = pl.read_csv("results/" + file, separator=",")
    axs[0].plot(df["eta"], df["y1"], label=f"Iter {file[-5]}")
    axs[1].plot(df["eta"], df["y2"], label=f"Iter {file[-5]}")
    axs[2].plot(df["eta"], df["y3"], label=f"Iter {file[-5]}")

    axs[0].set_ylabel(r"$f(\eta)$")
    axs[1].set_ylabel(r"$f'(\eta)$")
    axs[2].set_ylabel(r"$f''(\eta)$")

    axs[0].set_xlabel(r"$\eta$")
    axs[1].set_xlabel(r"$\eta$")
    axs[2].set_xlabel(r"$\eta$")

    axs[0].legend()
    axs[1].legend()
    axs[2].legend()

    fig.suptitle("Perfis de similaridade da solução de blasius", fontsize=24)


fig.savefig("perfis_de_similaridade.png", bbox_inches="tight", dpi=300)
