import argparse

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Anima um heatmap 2D projetado da barra 1D"
    )
    _ = parser.add_argument("--arquivo", type=str, required=True)
    _ = parser.add_argument("--tempo_total", type=float, required=True)
    _ = parser.add_argument("--dt", type=float, required=True)
    _ = parser.add_argument("--pontos", type=int, required=True)
    _ = parser.add_argument("--titulo", type=str, default="Heatmap da Barra")
    args = parser.parse_args()

    n = args.pontos
    dt = args.dt

    try:
        dados = np.fromfile(args.arquivo, dtype=np.float64)
    except FileNotFoundError:
        print(f"Erro: Arquivo {args.arquivo} não encontrado.")
        return

    n_frames = len(dados) // n
    dados = dados.reshape((n_frames, n))

    # Cria uma malha 2D a partir do vetor 1D
    # Queremos altura = 2 * x_star. x_star vai de -1 a 1, então a altura total é 2.
    # Vamos criar uma malha 2D onde cada coluna tem o mesmo valor da temperatura na posição x.
    height_points = 20  # Resolução vertical do heatmap
    heatmap_data = np.repeat(dados[:, np.newaxis, :], height_points, axis=1)

    fig, ax = plt.subplots(figsize=(12, 4))

    # limites: x de -1 a 1, y de -1 a 1
    im = ax.imshow(
        heatmap_data[0],
        extent=(-1, 1, -1, 1),
        origin="lower",
        aspect="auto",
        cmap="hot",
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Temperatura (°C)")

    ax.set_title(args.titulo)
    ax.set_xlabel("x* (Adimensional)")
    ax.set_ylabel("Altura")

    text_tempo = ax.text(
        0.02, 0.9, "", transform=ax.transAxes, color="Black", fontweight="bold"
    )

    def init():
        im.set_data(heatmap_data[0])
        text_tempo.set_text("")
        return [im, text_tempo]

    def update(frame):
        im.set_data(heatmap_data[frame])
        # Ajusta os limites de cor dinamicamente ou mantém fixo?
        # Para melhor visualização da evolução, vamos fixar no range global se possível
        text_tempo.set_text(f"Tempo: {frame * dt:.2f} s")
        return [im, text_tempo]

    # Ajusta os limites de cor para o range total de temperatura
    im.set_clim(vmin=np.min(dados), vmax=np.max(dados))

    # Para 10 segundos de vídeo:
    video_duration = 15.0  # s
    fps = 20
    total_frames_target = int(video_duration * fps)

    step = max(1, n_frames // total_frames_target)
    frames_to_show = range(0, n_frames, step)

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=frames_to_show,
        init_func=init,
        blit=True,
        interval=1000 / fps,
    )

    output_name = "heatmap_animado.mp4"
    print(f"Salvando heatmap como {output_name}...")
    try:
        ani.save(output_name, writer="ffmpeg", fps=fps, dpi=300)
    except:
        ani.save("heatmap_animado.gif", writer="pillow", fps=20)

    # plt.show()


if __name__ == "__main__":
    main()
