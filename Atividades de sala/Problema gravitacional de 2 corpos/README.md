# Problema Gravitacional de 2 Corpos

Esta atividade consiste na simulação numérica da órbita de dois corpos sob atração gravitacional mútua. A simulação compara dois métodos numéricos para resolver o sistema de equações diferenciais ordinárias (EDOs) de segunda ordem derivado da Lei da Gravitação Universal de Newton: o **Método de Euler** e o **Método Leapfrog**.

---

## Análise

Abaixo estão as respostas para as questões propostas sobre a conservação de energia e a estabilidade orbital dos dois métodos numéricos analisados.

### 1. O método de Euler conserva a energia total do sistema?

**Não.** O método de Euler não conserva a energia total do sistema. Durante a simulação, a energia total aumenta de forma acumulativa e monótona. No teste com $\Delta t = 1.0\times 10^{-3}$, a energia total do sistema iniciou em $-0.75$ e subiu para aproximadamente $-0.60$ ao fim de $t = 10$ (uma variação significativa de cerca de $20\%$).

### 2. O que acontece com as órbitas quando o método de Euler é usado por muitos passos de tempo?

As órbitas sofrem um efeito de espiralamento para fora (expansão orbital). Conforme o tempo passa e a energia artificial do sistema aumenta, os corpos ganham velocidade e se afastam cada vez mais de suas trajetórias elípticas fechadas originais, desintegrando a estabilidade física do sistema a longo prazo.

### 3. O método Leapfrog apresenta melhor comportamento na conservação de energia?

**Sim, muito superior.** O método Leapfrog preserva a energia mecânica de forma excelente. Com o mesmo passo de tempo $\Delta t = 1.0\times 10^{-3}$, a energia inicial de $-0.75$ terminou em aproximadamente $-0.749999984$ após 10.000 passos (um desvio residual irrelevante de apenas $\approx 1.6 \times 10^{-8}$). Além de extremamente pequeno, o erro de energia no método Leapfrog não diverge de forma monótona; ele oscila de forma estritamente limitada ao longo do tempo.

### 4. Qual método é mais adequado para esse problema conservativo? Justifique

O método **Leapfrog** é muito mais adequado. O problema gravitacional de dois corpos é um sistema mecânico conservativo (Hamiltoniano). O método Leapfrog é um **integrador simplético** de segunda ordem, o que significa que ele preserva a geometria de área no espaço de fase (teorema de Liouville). Essa propriedade simplética garante que a energia do sistema oscile estritamente em torno do valor real em vez de sofrer desvios sistemáticos (drift), mantendo a estabilidade orbital por períodos indefinidos de simulação. O método de Euler, por ser de primeira ordem e não simplético, é inadequado para integrações de longo prazo em sistemas físicos conservativos.

---

## 🛠️ Estrutura do Código

### 1. [main.cpp](./main.cpp)

O código em C++ está estruturado em classes reutilizáveis com as seguintes responsabilidades:

* `estado`: Estrutura simples para armazenar coordenadas $x$ e $y$ (posições, velocidades, acelerações, etc.).
* `configuracao_simulacao`: Gerencia os parâmetros do sistema físico ($G, m_1, m_2$), condições iniciais, passo de tempo ($\Delta t$), tempo total de simulação e o método de integração. Permite a configuração interativa via terminal.
* `integrador`: Contém a lógica física e matemática.
  * `calcular_derivadas_v()`: Calcula as acelerações baseando-se nas posições atuais através da força gravitacional mutua.
  * `dar_passo()`: Atualiza as posições e velocidades das partículas de acordo com o método selecionado (`euler` ou `leapfrog`).
* `solver`: Inicializa os vetores de estado temporais, gerencia o laço principal da simulação, calcula a energia mecânica total ($E = K + U$) em cada ponto, e salva os resultados finais em formato `.csv` no diretório `./results/`.

### 2. [plot.py](plot.py)

Script em Python que carrega os dados gerados pelos arquivos CSV de simulação, plota as trajetórias espaciais e o espaço de fase de ambos os corpos, salvando os gráficos comparativos:

* `Resultados_Euler.png`
* `Resultados_Leapfrog.png`

---

## Como Executar

### Pré-requisitos

* Compilador C++ compatível com C++17 (ex: `g++`)
* Python 3 com as bibliotecas `polars` e `matplotlib` instaladas.

### Passos para compilação e execução

1. **Compilar o programa em C++**:

    ```bash
    g++ -O3 main.cpp -o main
    ```

2. **Executar as simulações**:

    ```bash
    ./main
    ```

    Isso gerará os arquivos `resultados_euler.csv` e `resultados_leapfrog.csv` dentro da pasta `./results/`.

3. **Gerar os gráficos comparativos**:

    ```bash
    python plot.py
    ```

    Isso salvará as imagens `Resultados_Euler.png` e `Resultados_Leapfrog.png` no diretório raiz.
