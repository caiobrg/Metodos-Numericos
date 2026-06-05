import argparse

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description="Anima a validacao final")
    _ = parser.add_argument(
        "--arquivo", type=str, required=True, help="Arquivo de resultados (results.bin)"
    )
    _ = parser.add_argument("--tempo_total", type=float, required=True)
    _ = parser.add_argument("--dt", type=float, required=True)
    _ = parser.add_argument("--pontos", type=int, required=True)
    _ = parser.add_argument("--L", type=float, required=True)
    _ = parser.add_argument("--k", type=float, required=True)
    _ = parser.add_argument("--cp", type=float, required=True)
    _ = parser.add_argument("--rho", type=float, required=True)
    _ = parser.add_argument("--h1", type=float, required=True)
    _ = parser.add_argument("--t_inf", type=float, required=True)
    _ = parser.add_argument("--t_cond_inicial", type=float, required=True)
    args = parser.parse_args()

    n = args.pontos
    dt = args.dt
    L = args.L
    k = args.k
    cp = args.cp
    rho = args.rho
    h1 = args.h1
    t_inf = args.t_inf
    t_cond_inicial = args.t_cond_inicial

    # Parâmetros físicos para a solução analítica
    alpha = k / (cp * rho)
    L_meio = L / 2.0
    Bi_ana = (h1 * L_meio) / k

    # Lê o arquivo binário
    try:
        dados_completos = np.fromfile(args.arquivo, dtype=np.float64)
    except FileNotFoundError:
        print(f"Erro: Arquivo {args.arquivo} não encontrado.")
        return

    # O número de frames é o total de elementos / n
    n_frames = len(dados_completos) // n
    dados_completos = dados_completos.reshape((n_frames, n))

    # Cálculo dos autovalores
    def f_lambda(lb, Bi):
        """
        Função que devemos achar os 0s
        """
        return lb * np.tan(lb) - Bi

    def calcular_autovalores(Bi, num_termos=100):
        """
        Calcula os autovalores pelo método da bissecção
        """
        lambdas = []
        for i in range(1, num_termos + 1):
            a = (i - 1) * np.pi + 1e-6
            b = (i - 1) * np.pi + np.pi / 2.0 - 1e-6
            # Bisseção simples
            for _ in range(50):
                c = (a + b) / 2.0
                if f_lambda(c, Bi) * f_lambda(a, Bi) < 0:
                    b = c
                else:
                    a = c
            lambdas.append((a + b) / 2.0)
        return np.array(lambdas)

    lambdas = calcular_autovalores(Bi_ana, 500)

    def theta_analitico(x_star, Fo, lbs):
        """
        Somatório para a temperatura admensional
        """
        if Fo <= 0:
            return 1.0
        coef = (4.0 * np.sin(lbs)) / (2.0 * lbs + np.sin(2.0 * lbs))
        termos = coef * np.cos(lbs * x_star) * np.exp(-(lbs**2) * Fo)
        return np.sum(termos)

    # Malha x_star
    dx = L / (n - 1)
    x_star = np.array([((i * dx) - L_meio) / L_meio for i in range(n)])

    # Prepara a figura
    fig, ax = plt.subplots(figsize=(10, 6))
    (line_num,) = ax.plot([], [], label="Numérica (Diferenças Finitas)", lw=2)
    (line_ana,) = ax.plot([], [], linestyle="--", label="Analítica (Fourier)", lw=2)
    text_r2 = ax.text(
        0.05,
        0.9,
        "",
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        bbox=dict(facecolor="white", alpha=0.5),
    )
    text_tempo = ax.text(0.05, 0.85, "", transform=ax.transAxes, fontsize=12)

    ax.set_xlim(-1.1, 1.1)
    # Ajusta o limite Y com base nos dados iniciais e finais
    y_min = min(t_inf, t_cond_inicial) - 10
    y_max = max(t_inf, t_cond_inicial) + 10
    ax.set_ylim(y_min, y_max)

    ax.set_title(
        "Validação Transiente - Evolução Temporal", fontsize=14, fontweight="bold"
    )
    ax.set_xlabel("Posição Adimensional ($x^*$)", fontsize=12)
    ax.set_ylabel(r"Temperatura ($^\circ$C)", fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    def init():
        line_num.set_data([], [])
        line_ana.set_data([], [])
        text_r2.set_text("")
        text_tempo.set_text("")
        return line_num, line_ana, text_r2, text_tempo

    def update(frame):
        tempo = frame * dt
        T_num = dados_completos[frame]

        Fo = (alpha * tempo) / (L_meio**2)
        T_ana = np.array(
            [
                t_inf + theta_analitico(xs, Fo, lambdas) * (t_cond_inicial - t_inf)
                for xs in x_star
            ]
        )

        # R^2
        media_ana = np.mean(T_ana)
        ss_res = np.sum((T_num - T_ana) ** 2)
        ss_tot = np.sum((T_ana - media_ana) ** 2)
        r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

        line_num.set_data(x_star, T_num)
        line_ana.set_data(x_star, T_ana)
        text_r2.set_text(f"$R^2 = {r2:.6f}$")
        text_tempo.set_text(f"Tempo = {tempo:.2f} s")

        return line_num, line_ana, text_r2, text_tempo

    # Para 10 segundos de vídeo:
    video_duration = 15.0  # s
    fps = 20
    total_frames_target = int(video_duration * fps)

    # Pula alguns frames para a animação não ser muito lenta se houver muitos passos
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

    # Salva a animação
    print("Salvando animação como 'validacao_animada.mp4'...")
    try:
        ani.save("validacao_animada.mp4", writer="ffmpeg", fps=fps, dpi=150)
        print("Sucesso! Animação salva.")
    except Exception as e:
        print(f"Erro ao salvar mp4: {e}. Tentando salvar como GIF...")
        try:
            ani.save("validacao_animada.gif", writer="pillow", fps=20)
            print("Sucesso! Animação salva como GIF.")
        except Exception as e2:
            print(f"Erro ao salvar GIF: {e2}")

    # plt.show()


if __name__ == "__main__":
    main()
