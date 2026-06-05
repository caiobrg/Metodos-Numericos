# PPC3 - Sistemas Lineares

## Resumo

Este projeto busca implementar um solver numérico para resolver a equação do calor unidimensional transiente utilizando o método de diferenças finitas e o Algoritmo de Thomas para sistemas de equações lineares tridiagonais. O projeto foi estruturado da seguinte maneira:

``` plain text
└── 📁PPC3 - Sistemas lineares
    ├── .gitignore
    ├── heatmap_animation.py
    ├── main.cpp
    ├── plot_convergencia.py
    ├── plot_validacao.py
    ├── README.md
    └── solver_diferencas_finitas.hpp
```

O motor do solver é o arquivo [solver_diferencas_finitas.hpp](solver_diferencas_finitas.hpp), que é um arquivo de cabeçalho contendo a implementação do método de diferenças finitas e a resolução do sistema linear gerado pelo algoritmo de Thomas. O arquivo orquestrador utilizado para realizar as simulações e análises é o [main.cpp](main.cpp). Por fim, para realizar a plotagem da convergência e as animações do perfil de temperatura e do mapa de calor, elaboraram-se scripts em Python: [plot_convergencia.py](plot_convergencia.py), [plot_validacao.py](plot_validacao.py) e [heatmap_animation.py](heatmap_animation.py).

Caso deseje somente as instruções de compilação e uso, favor se dirigir à seção [Instruções de uso](#instruções-de-uso).

## 1 Introdução

A análise da transferência de calor transiente em uma dimensão é fundamental em diversas aplicações de engenharia, como no estudo do resfriamento e aquecimento de componentes. A equação diferencial parcial que governa a condução de calor unidimensional com geração interna de energia térmica é dada por:

$$\frac{\partial^2 T}{\partial x^2} + \frac{\dot{q}}{k} = \frac{1}{\alpha} \frac{\partial T}{\partial t}$$

Onde $k$ é a condutividade térmica, $\alpha$ é a difusividade térmica e $\dot{q}$ é a taxa de geração interna de calor.

Para resolver este problema computacionalmente, podemos discretizar o domínio no espaço e no tempo através do método de diferenças finitas. O processo gera um sistema de equações algébricas lineares a ser resolvido em cada instante de tempo. Como as interações térmicas ocorrem apenas entre nós adjacentes (vizinhos imediatos), a matriz de coeficientes resultante é tridiagonal.

Aqui está uma versão aprimorada da **Seção 2**, com uma explicação mais detalhada sobre como o método das diferenças finitas é aplicado na aproximação das derivadas e qual o impacto real da introdução da ponderação $\beta$. A Seção 3 foi mantida com o texto que você validou.

## 2 Implementação do solver

O método das Diferenças Finitas consiste em discretizar o domínio contínuo do problema em uma malha com espaçamento espacial $\Delta x$ e calcular a evolução térmica em passos discretos de tempo $\Delta t$. Para isso, as derivadas da equação do calor são substituídas por aproximações algébricas: a derivada temporal (primeira ordem) é aproximada por uma diferença avançada, enquanto a derivada espacial (difusão, de segunda ordem) é aproximada por diferenças centrais relacionando a temperatura de um nó $i$ à dos seus vizinhos imediatos $i-1$ e $i+1$.

**O papel do parâmetro $\beta$**
Na formulação discreta para um nó interno genérico, surge a questão: os termos espaciais de difusão de calor devem ser avaliados no tempo atual (conhecido) ou no tempo futuro (desconhecido)? O parâmetro $\beta \in [0,1]$ é introduzido justamente para ponderar essas duas avaliações, generalizando o método:

* **$\beta = 0$ (Euler Explícito):** A difusão é avaliada inteiramente no tempo atual. As equações são fáceis de resolver individualmente, mas o método é condicionalmente estável, exigindo passos de tempo ($\Delta t$) muito pequenos para que a solução não divirja.
* **$\beta = 1$ (Implícito Puro):** A difusão é avaliada no tempo futuro. O método se torna incondicionalmente estável, o que permite adotar saltos de tempo maiores, ao custo de exigir a resolução de um sistema de equações lineares a cada passo.
* **$\beta = 0.5$ (Crank-Nicolson):** Faz uma média aritmética entre os instantes atual e futuro. É um esquema muito popular por equilibrar a estabilidade incondicional do implícito puro e a precisão do método de Euler Explícito.

Ao agrupar todos os termos avaliados no instante futuro ($t + \Delta t$) no lado esquerdo da equação e os termos no instante atual ($t$) no lado direito, o sistema de equações algébricas assume a forma matricial:

$$\underbar{\underbar{A}} \cdot \underbar{T} = \underbar{B}$$

sendo $\underbar{\underbar{A}}$ uma matriz tridiagonal, uma vez que cada nó da malha interage apenas com seus vizinhos adjacentes.

A resolução de sistemas tridiagonais utilizando métodos diretos convencionais (como a Eliminação de Gauss genérica) é custosa e computacionalmente ineficiente. Por isso, a implementação em `solver_diferencas_finitas.hpp` utiliza o **Algoritmo de Thomas** (Tridiagonal Matrix Algorithm - TDMA). O método otimiza a eliminação Gaussiana para matrizes tridiagonais com um custo operacional drasticamente reduzido a $O(N)$ operações por passo temporal, garantindo alta eficiência por meio de uma etapa de decomposição progressiva e outra de substituição regressiva.

## 3 Estudo de Convergência

O orquestrador em C++ implementa rotinas para estudar a convergência numérica frente à discretização da malha espacial (variando o número de nós $N$ para um $\Delta t$ fixo) e temporal (variando o $\Delta t$ para um $N$ fixo).

A verificação é feita confrontando a solução numérica do transiente de temperatura da placa (sem geração de calor) com a solução analítica exata, obtida via série de Fourier e raízes transcendentais baseadas no Número de Biot ($Bi$). O R-quadrado ($R^2$) ou o erro $(1 - R^2)$ são utilizados como métricas. O script `plot_convergencia.py` é responsável por ler os resultados numéricos e exibir os gráficos de convergência espacial e temporal em escala log-log, indicando a ordem do erro.

## 4 Validação Animada e Heatmap

Para verificar fisicamente os resultados, o programa realiza simulações animadas do comportamento transiente da temperatura.

1. **Validação Animada:** O código realiza a comparação dinâmica da evolução do perfil de temperaturas $T(x,t)$ na placa entre a solução numérica e a solução analítica (via série de Fourier) até o equilíbrio térmico. Um vídeo de validação é exportado através do `plot_validacao.py`.
2. **Heatmap:** O orquestrador também implementa uma análise de caso mais complexa que inclui a inserção de um termo de geração interna volumétrica na placa, representando um processo físico de dissipação ou reação. O transiente é visualizado sob a forma de um mapa de calor (Heatmap) que se altera no tempo, processado em `heatmap_animation.py`.

## Instruções de Uso

O projeto utiliza um motor em C++ e scripts de visualização em Python. Abaixo estão as etapas para gerar os resultados.

### Pré-requisitos

* **Compilador C++**: GCC (`g++`) ou Clang com suporte à C++11 ou superior.
* **Python 3.x**: Recomendado instalação via `pip` das bibliotecas `numpy`, `pandas`, `matplotlib`.
* **FFmpeg**: Requerido pelo matplotlib para exportar as animações transientes para `.mp4`. Se ausente, o código buscará exportar como `.gif`.

### Passo 1: Compilação do Orquestrador (C++)

O orquestrador numérico (`main.cpp`) deve estar no mesmo diretório do arquivo `solver_diferencas_finitas.hpp`. Compile o código ativando otimizações:

```bash
g++ -O3 main.cpp -o main
```

### Passo 2: Execução

Ao executar o orquestrador, ele cuidará de criar as pastas de resultados, chamar as funções do solver, gerar arquivos em disco (`.bin` e `.csv`) e invocar automaticamente os scripts Python associados utilizando `std::system`.

```bash
./main
```

Durante a execução, você verá na tela a evolução dos erros computados das rotinas de convergência, e a notificação da geração das figuras e das animações.

* **Configurações:** Você pode alterar propriedades físicas e numéricas do material na função `main()` como $L$ (comprimento), $dt$ (passo de tempo), $\beta$ (0.5 para Crank-Nicolson), as propriedades térmicas ($k, c_p, \rho$) e a malha de elementos diretamente no `main.cpp` e recompilar.

### Passo 3: Verificação de Saídas

Os resultados poderão ser conferidos nos seguintes arquivos:

1. `estudo_convergencia.png`: Imagem contendo as curvas de convergência espacial e temporal indicando a ordem calculada.
2. `validacao_animada.mp4` (ou `.gif`): Vídeo comparando a curva de difusão de calor do solver ao da solução analítica no tempo.
3. `heatmap_animado.mp4` (ou `.gif`): Vídeo do mapa de cor da temperatura da placa (simulando difusão + geração de calor).
4. `results/`: Diretório contendo as extrações temporárias em `bin` e os dados agregados da convergência espacial e temporal em arquivos `csv`.
