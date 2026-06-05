import polars as pl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


def plot_fractal_shaded(filename, output_name, cmap="jet"):
    print(f"Lendo e renderizando {filename}...")
    try:
        # df = pd.read_csv(filename)
        df = pl.read_csv(filename)
    except FileNotFoundError:
        print(f"Arquivo {filename} não encontrado.")
        return

    if not df.filter(pl.col("sucesso") == 1).is_empty():
        # Cria colunas temporárias arredondadas apenas para achar o ID da bacia
        unique_basins = (
            df.filter(pl.col("sucesso") == 1)
            .select(
                [
                    pl.col("r_conv").round(2).alias("r_round"),
                    pl.col("s_conv").round(2).alias("s_round"),
                ]
            )
            .unique()
            .with_row_index(name="basin_id")
        )

        # Na hora de juntar, arredonda temporariamente as colunas do df principal no próprio left_on
        df = df.join(
            unique_basins,
            left_on=[pl.col("r_conv").round(2), pl.col("s_conv").round(2)],
            right_on=["r_round", "s_round"],
            how="left",
        )

        # 4. Tratamos os casos que não convergiram (nulos) colocando -1
        df = df.with_columns(
            pl.col("basin_id").fill_null(-1).alias("basin")
        ).drop("basin_id")
    else:
        df = df.with_columns(pl.col().fill_null(-1).alias("basin"))

    # 1. Faz o pivot e ORDENA o eixo Y (linhas "s")
    pivot_basin = df.pivot(index="s", on="r", values="basin").sort("s")
    pivot_iter = df.pivot(index="s", on="r", values="iteracoes").sort("s")
    pivot_success = df.pivot(index="s", on="r", values="sucesso").sort("s")
    pivot_exact = df.pivot(index="s", on="r", values="exact_root").sort("s")

    # 2. Extrai os nomes das colunas (que o Polars transformou em string),
    # ignora a coluna "s", e ordena numericamente (convertendo para float na chave de ordenação)
    r_cols = sorted(
        [col for col in pivot_basin.columns if col != "s"], key=float
    )

    # Recria a lista de colunas na ordem perfeita (eixo X da esquerda pra direita)
    ordered_cols = ["s"] + r_cols

    # 3. Reorganiza as colunas dos DataFrames
    pivot_basin = pivot_basin.select(ordered_cols)
    pivot_iter = pivot_iter.select(ordered_cols)
    pivot_success = pivot_success.select(ordered_cols)
    pivot_exact = pivot_exact.select(ordered_cols)

    # 4. Agora sim, extrai para matrizes NumPy removendo o "s"
    basin_grid = pivot_basin.drop("s").to_numpy()
    iter_grid = pivot_iter.drop("s").to_numpy()
    success_grid = pivot_success.drop("s").to_numpy()
    exact_grid = pivot_exact.drop("s").to_numpy()
    # Configuração do Colormap
    num_basins = int(df["basin"].max() + 1)
    # 1. Pegamos a paleta original completa (ex: 'magma' ou 'inferno')
    original_cmap = plt.colormaps.get_cmap(cmap)

    # 2. Criamos um vetor de 0.3 a 1.0 (Pulando os primeiros 30% da paleta)
    # Se ainda estiver escuro, aumente o 0.3 para 0.4. Se ficar muito claro, baixe para 0.2.
    color_range = np.linspace(0.3, 1.0, 256)

    # 3. Extraímos apenas as cores desse intervalo e criamos um novo colormap
    colors_subset = original_cmap(color_range)
    cmap = mcolors.ListedColormap(colors_subset)

    # 1. Componente de cor base em RGB (em vez de HSV)
    rgb_image = np.zeros((*basin_grid.shape, 3))

    if num_basins > 0:
        for i in range(num_basins):
            # Normaliza o ID da bacia entre 0.0 e 1.0 para pegar a cor correta no gradiente
            norm_val = i / (num_basins - 1) if num_basins > 1 else 0.5
            # O cmap retorna (R, G, B, Alpha), pegamos apenas os 3 primeiros [:3]
            rgb_image[basin_grid == i] = cmap(norm_val)[:3]

    # 2. Máscara de Sombreamento (Shading)
    max_iter = iter_grid.max() if iter_grid.max() > 0 else 1
    dy, dx = np.gradient(iter_grid)
    light = np.array([-1, 1, 0.5])
    light = light / np.linalg.norm(light)

    n = np.dstack(
        [
            -dx * 100.0 / max_iter,
            -dy * 100.0 / max_iter,
            np.ones_like(iter_grid),
        ]
    )
    n = n / np.linalg.norm(n, axis=2, keepdims=True)

    shading = np.sum(n * light, axis=2)
    shading = np.clip(shading, 0, 1)

    # 3. Mistura a cor RGB com o sombreamento e o decaimento
    shading = shading * 0.5 + 0.5
    iter_decay = 1.0 - (iter_grid / max_iter) ** 0.3

    # Calcula o multiplicador final de brilho/sombra
    brightness_mask = shading * (iter_decay * 0.4 + 0.6) * 1.1

    # Multiplica cada canal de cor (R, G e B) pela máscara de brilho
    rgb_image = rgb_image * np.expand_dims(brightness_mask, axis=2)

    # Garante que as cores não passem do limite de exibição [0, 1]
    rgb_image = np.clip(rgb_image, 0, 1)

    # 4. Pós-processamento de pontos especiais
    rgb_image[success_grid == 0] = [0, 0, 0]
    rgb_image[exact_grid == 1] = [1, 1, 1]

    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(
        rgb_image,
        extent=[df["r"].min(), df["r"].max(), df["s"].min(), df["s"].max()],
        origin="lower",
    )

    # Adicionando legenda para as bacias
    if num_basins > 0:
        from matplotlib.lines import Line2D

        legend_elements = []
        # Encontrar um representante para cada bacia rapidamente
        representatives = (
            df.filter(
                pl.col("basin") >= 0
            )  # Filtra apenas os que convergiram (ignora os -1)
            .select(
                ["r_conv", "s_conv", "basin"]
            )  # Pega só as colunas que importam
            .unique(subset=["basin"])  # Não seleciona valores duplicados
            .sort("basin")  # Ordenamento
        )

        for rep in representatives.iter_rows(named=True):
            i = rep["basin"]

            # Normaliza o i exatamente como fizemos na imagem
            norm_val = i / (num_basins - 1) if num_basins > 1 else 0.5
            color = cmap(norm_val)[:3]

            label = (
                f"Par de Raízes: r={rep['r_conv']:.2f}, s={rep['s_conv']:.2f}"
            )

            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    label=label,
                    markerfacecolor=color,
                    markersize=10,
                )
            )

        # Adiciona entrada para não convergente
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Não Convergiu",
                markerfacecolor="black",
                markersize=10,
            )
        )

        # Adiciona entrada para convergência imediata
        legend_elements.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label="Conv. Imediata",
                markerfacecolor="white",
                markersize=10,
            )
        )

        ax.legend(
            handles=legend_elements,
            loc="upper right",
            bbox_to_anchor=(1.45, 1.0),
        )

    ax.set_title(f"Bairstow Fractal - {filename}")
    ax.set_xlabel("r")
    ax.set_ylabel("s")
    plt.savefig(
        output_name,
        dpi=400,
        bbox_inches="tight",
        pil_kwargs={
            "optimize": True,
            "quality": 80,
        },  # O 'quality' afeta JPG e WebP
    )
    plt.close()


plot_fractal_shaded(
    "outputs/fractal_pol_carac.csv", "figures/fractal_pol_carac.webp"
)
plot_fractal_shaded("outputs/fractal_pol_7.csv", "figures/fractal_pol_7.webp")
plot_fractal_shaded(
    "outputs/fractal_pol_comp.csv", "figures/fractal_pol_comp.webp"
)
