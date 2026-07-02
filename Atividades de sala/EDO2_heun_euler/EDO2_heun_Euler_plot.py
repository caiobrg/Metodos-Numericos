import matplotlib.pyplot as plt
import polars as pl

df_euler = pl.read_csv("integração_euler.dat", separator=",")
df_heun = pl.read_csv("integração_heun.dat", separator=",")

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
ax[0].plot(df_euler["t"], df_euler["x"], label="Euler")
ax[0].plot(df_heun["t"], df_heun["x"], label="Heun")
ax[0].set_xlabel("Tempo [s]")
ax[0].set_ylabel("Posição [m]")
ax[0].legend()

ax[1].plot(df_euler["x"], df_euler["v"], label="Euler")
ax[1].plot(df_heun["x"], df_heun["v"], label="Heun")
ax[1].set_xlabel("Posição [m]")
ax[1].set_ylabel("Velocidade [m/s]")
ax[1].set_aspect(abs(df_euler["x"].to_numpy().max() / df_euler["v"].to_numpy().max()))
ax[1].legend()

fig.savefig("integração_euler.png", bbox_inches="tight", dpi=300)
