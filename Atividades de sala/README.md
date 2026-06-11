# Atividades de Sala - Métodos Numéricos

Este diretório contém implementações de algoritmos e métodos numéricos desenvolvidos em C++ para resolver problemas de otimização, localização de raízes e resolução de sistemas lineares. Abaixo está uma breve explicação de cada código, com foco em sua funcionalidade geral e variáveis-chave.

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
