import matplotlib.pyplot as plt
import polars as pl


def plot(df: pl.DataFrame, filename: str) -> None:
    fig, axs = plt.subplots(ncols=2, figsize=(18, 6))
    axs[0].plot(df["x1"], df["y1"], label="Corpo 1")
    axs[0].plot(df["x2"], df["y2"], label="Corpo 2")

    axs[0].set_title("Trajetórias")
    axs[0].legend()

    axs[1].plot(
        df["x1"] / df["x1"].max(),
        df["vx1"] / df["vx1"].max(),
        label="Corpo 1",
    )
    axs[1].plot(
        df["x2"] / df["x2"].max(),
        df["vx2"] / df["vx2"].max(),
        label="Corpo 2",
    )
    axs[1].set_title("Espaço de fase")
    axs[1].legend()
    fig.savefig(filename, bbox_inches="tight", dpi=300)


df_euler = pl.read_csv("./results/resultados_euler.csv")
df_leap = pl.read_csv("./results/resultados_leapfrog.csv")

plot(df_euler, "Resultados_Euler.png")
plot(df_leap, "Resultados_Leapfrog.png")
