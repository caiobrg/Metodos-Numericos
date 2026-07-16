# Atividades de Sala - Métodos Numéricos

Este diretório contém implementações de algoritmos e métodos numéricos desenvolvidos em C++ para resolver problemas de otimização, localização de raízes, resolução de sistemas lineares, equações diferenciais ordinárias (EDOs) e equações diferenciais parciais (EDPs). Abaixo está uma breve explicação de cada código, com foco em sua funcionalidade geral e variáveis-chave.

---

## 📌 Sumário dos Códigos

### 1. `aurea.cpp` - Método da Seção Áurea

* **Funcionalidade Geral:** Realiza otimização unidimensional para encontrar o valor de resistência ajustável ($R_a$) que maximiza a potência dissipada ($P$) em um circuito elétrico. Utiliza o método de busca da seção áurea, estreitando iterativamente o intervalo contendo o ponto máximo.

* **Variáveis-Chave:**
  * `ra_lower` e `ra_upper`: Limites inferior e superior do intervalo de busca para $R_a$.
  * `razao`: A constante da razão áurea ($\frac{\sqrt{5}-1}{2} \approx 0.618033$) usada para particionar o intervalo.
  * `ra1` e `ra2`: Pontos internos de amostragem gerados pela razão áurea.
  * `ra_opt` e `p_max`: Resistência ótima encontrada e potência máxima calculada.
  * `tol` e `erro`: Tolerância de parada e erro estimado.

### 2. `decomp_LU.cpp` - Decomposição LU Manual

* **Funcionalidade Geral:** Resolve um sistema de equações lineares de ordem 10 ($Ax = B$) gerado de forma aleatória implementando o algoritmo de Decomposição LU e as substituições (progressiva e regressiva) manualmente através de laços de repetição `for`, sem utilizar resolvedores embutidos da biblioteca `Eigen`. Inclui a validação dos resultados obtidos por meio do cálculo das normas dos erros residuais e da verificação direta da igualdade $A = LU$.
* **Variáveis-Chave:**
  * `A` e `B`: Matriz de coeficientes aleatória ($10 \times 10$) e vetor de termos independentes aleatórios, respectivamente.
  * `L` e `U`: Matrizes fatoradas triangular inferior (Lower, com diagonal unitária) e superior (Upper) calculadas manualmente.
  * `Y`: Vetor intermediário obtido na substituição progressiva ($LY = B$).
  * `X`: Vetor solução final obtido na substituição regressiva ($UX = Y$).

### 3. `interpolacao.cpp` - Método de Interpolação Quadrática

* **Funcionalidade Geral:** Otimização unidimensional utilizando interpolação polinomial quadrática (de 3 pontos) para determinar o máximo da função de potência. A cada iteração, ajusta-se uma parábola aos três pontos atuais, calcula-se o vértice para encontrar uma nova estimativa e atualizam-se os limites dinamicamente conforme os valores de potência obtidos.
* **Variáveis-Chave:**
  * `ra_0`, `ra_1` e `ra_2`: Os três valores de estimativa de resistência mantidos e atualizados a cada passo.
  * `ra_3`: Nova estimativa calculada no vértice da parábola interpoladora.
  * `p_0`, `p_1`, `p_2` e `p_3`: Potências calculadas nos respectivos pontos de resistência.
  * `calculate_ra3`: Função lambda que implementa a fórmula matemática para encontrar a abscissa do vértice da parábola.
  * `erro` e `tol`: Critério de parada baseado na variação relativa entre estimativas sucessivas de $R_a$.

### 4. `muller.cpp` - Método de Muller

* **Funcionalidade Geral:** Localiza raízes de funções polinomiais não lineares de grau 3. O método é uma generalização do método da secante que traça uma parábola através de três pontos aproximados para projetar o cruzamento com o eixo horizontal. O resultado da raiz encontrada é salvo no terminal e em um arquivo de texto (`raizes.dat`).
* **Variáveis-Chave:**
  * `x0`, `x1` e `x2`: Três chutes iniciais para o algoritmo.
  * `x3`: A estimativa atualizada da raiz obtida a partir da parábola de interpolação.
  * `h0`, `h1`, `delta0` e `delta1`: Variáveis de passo e razões de diferenças divididas.
  * `a`, `b` e `c`: Coeficientes do polinômio quadrático de ajuste local ($ay^2 + by + c = 0$).
  * `discriminante`: Valor sob o radical para determinação do passo com maior proximidade da raiz.

### 5. `newton_raphson.cpp` - Método de Newton-Raphson

* **Funcionalidade Geral:** Encontra as raízes de uma equação algébrica não linear usando o método de aproximações sucessivas de Newton-Raphson clássico. Possui também uma função candidata para o Newton-Raphson Modificado (`f_NRmod`) voltado para tratamento de raízes múltiplas. Salva os resultados no arquivo `raizes.dat`.
* **Variáveis-Chave:**
  * `x_atual`: Estimativa de ponto flutuante da raiz ao longo das iterações.
  * `f(x)`, `f_derivada1(x)` e `f_derivada2(x)`: Computam, respectivamente, o polinômio, sua primeira derivada e segunda derivada.
  * `f_NR` e `f_NRmod`: Funções que calculam a próxima estimativa com base no passo clássico ou modificado.
  * `erro` e `tol`: Tolerância absoluta e diferença absoluta consecutiva usada para encerramento.

### 6. `EDO_metodo_do_tiro.cpp` - Método do Tiro para EDOs de Valor de Contorno

* **Funcionalidade Geral:** Resolve um problema de valor de contorno (BVP) para a equação de condução de calor unidimensional em uma barra, onde são conhecidas as temperaturas nas duas extremidades ($T(0)$ e $T(L)$). O método transforma o BVP em um problema de valor inicial (IVP), integrando via Método de Euler com diferentes chutes para a derivada inicial ($\frac{dT}{dx}\big|_{x=0}$). A cada iteração, o chute é refinado pelo **método da secante** até que a condição de contorno na extremidade $x = L$ seja satisfeita dentro da tolerância.
* **Variáveis-Chave:**
  * `L`, `N`, `q_dot`, `k`: Comprimento da barra [m], número de pontos da malha, geração interna de calor [W/m³] e condutividade térmica [W/(m·K)].
  * `T_0` e `T_L`: Condições de contorno (temperatura em $x = 0$ e $x = L$).
  * `f_1_ini` e `f_2_ini`: Dois chutes iniciais para a taxa de variação da temperatura ($\frac{dT}{dx}$).
  * `T_1` e `T_2`: Vetores de temperatura resultantes da integração com cada chute.
  * `f_3`: Novo chute calculado pela fórmula da secante com base nos erros de contorno.
  * `err1` e `err2`: Erros das condições de contorno ($T(L)_{\text{calculado}} - T_L$) para cada chute.

### 7. `EDO2_heun_euler/` - Integração de EDO de 2ª Ordem (Euler e Heun)

* **Funcionalidade Geral:** Simula um sistema massa-mola-amortecedor com forçamento harmônico ($m\ddot{x} + c\dot{x} + kx = F_0 \sin(\omega t)$) utilizando dois métodos de integração numérica para EDOs de segunda ordem: o **Método de Euler** e o **Método de Heun** (Euler melhorado / preditor-corretor). Os resultados são salvos em arquivos `.dat` para visualização com o script Python auxiliar.
* **Variáveis-Chave:**
  * `m`, `c`, `k`: Massa [kg], amortecimento viscoso [Ns/m] e rigidez [N/m] do sistema.
  * `f_0` e `w`: Amplitude [N] e frequência angular [rad/s] do forçamento harmônico.
  * `delta_t` e `t_total`: Passo de integração [s] e tempo total de simulação [s].
  * `x` e `v`: Vetores de deslocamento e velocidade ao longo do tempo.
  * `x_0` e `v_0` (Heun): Vetores de predição (Euler) usados como estimativa intermediária no passo corretor.
  * `dv_dt()`: Função que calcula a aceleração a partir da equação do movimento.
* **Arquivos:**
  * `EDO2_heun_Euler.cpp`: Código-fonte C++ com ambas as integrações.
  * `EDO2_heun_Euler_plot.py`: Script Python para plotagem dos resultados.
  * `integração_euler.png`: Gráfico gerado com os resultados da simulação.

### 8. `Problema gravitacional de 2 corpos/` - Simulação Gravitacional de 2 Corpos

* **Funcionalidade Geral:** Simula a órbita de dois corpos sob atração gravitacional mútua, comparando os métodos de **Euler** e **Leapfrog** para resolver o sistema de EDOs de segunda ordem derivado da Lei da Gravitação Universal de Newton. O código é estruturado em classes reutilizáveis com arquitetura orientada a objetos. Os resultados incluem posições, velocidades e energia mecânica total ($E = K + U$) ao longo do tempo, salvos em formato CSV.
* **Variáveis-Chave:**
  * `estado`: Estrutura contendo coordenadas $(x, y)$ para posições, velocidades e acelerações.
  * `configuracao_simulacao`: Classe que gerencia os parâmetros físicos ($G, m_1, m_2$), condições iniciais, passo de tempo ($\Delta t$) e método de integração. Configurável via terminal.
  * `integrador`: Classe que implementa a lógica de integração, com `calcular_derivadas_v()` (acelerações gravitacionais) e `dar_passo()` (atualização de posições e velocidades).
  * `solver`: Classe que gerencia o laço de simulação, calcula a energia mecânica total e exporta resultados para `.csv`.
* **Arquivos:**
  * `main.cpp`: Código-fonte C++ principal.
  * `plot.py`: Script Python para plotagem das trajetórias e espaço de fase.
  * `README.md`: Análise detalhada da conservação de energia e estabilidade orbital.
  * `Resultados_Euler.png` e `Resultados_Leapfrog.png`: Gráficos comparativos dos métodos.

### 9. `Distribuição de calor em placa/` - Distribuição de Temperatura em Placa 2D (Laplace)

* **Funcionalidade Geral:** Resolve a equação de Laplace bidimensional ($\nabla^2 T = 0$) para determinar a distribuição estacionária de temperatura em uma placa retangular (aleta) com condições de contorno de Dirichlet: temperatura $T_b$ na base e $T_{\infty}$ nas demais paredes. Utiliza o **método iterativo de Gauss-Seidel** aplicado ao esquema de diferenças finitas centradas de segunda ordem. O código é estruturado em classes orientadas a objetos com arquitetura modular e mede o tempo de execução do solver.
* **Variáveis-Chave:**
  * `configSimulacao`: Classe de configuração com parâmetros geométricos ($L, H$), malha ($N_x, N_y$), temperaturas de contorno ($T_b, T_\infty$), tolerância e número máximo de iterações.
  * `Node`: Classe representando um nó da malha, com identificador, posição $(x, y)$ e temperatura $T$.
  * `Mesh`: Classe responsável pela construção da malha de nós, aplicação das condições de contorno iniciais, mapeamento 2D→1D e exportação dos resultados em CSV.
  * `sistema`: Classe que implementa o solver iterativo de Gauss-Seidel, calculando em cada passo a média ponderada dos vizinhos (N, S, E, W) com fator $\beta = (\Delta x / \Delta y)^2$.
  * `beta`: Razão de aspecto da malha ($\beta = (\Delta x / \Delta y)^2$), usada nos coeficientes do esquema de diferenças finitas.
* **Arquivos:**
  * `main.cpp`: Código-fonte C++ principal com as classes e simulações para malhas de diferentes resoluções (11×11 até 500×500).
  * `plot.py`: Script Python para visualização dos campos de temperatura.
  * `results/`: Diretório com os arquivos CSV de malhas resolvidas.
