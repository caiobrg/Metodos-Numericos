# Métodos Numéricos - Cálculo Numérico Aplicado (CNA)

![Logo UnB](logo_unb.png)

[![C++](https://img.shields.io/badge/C%2B%2B-17%2F20-blue.svg?logo=c%2B%2B)](https://isocpp.org/)
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue?logo=python)](https://www.python.org/)
[![UnB](https://img.shields.io/badge/UnB-ENM-darkgreen)](https://enm.unb.br/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Este repositório contém os códigos, simulações e relatórios desenvolvidos ao longo da disciplina **Cálculo Numérico Aplicado (CNA)**, ofertada pelo **Departamento de Engenharia Mecânica (ENM)** da **Universidade de Brasília (UnB)** e ministrada pelo **Prof. Dr. Rafael Gabler Gontijo**.

Os projetos combinam motores numéricos implementados em **C++** com scripts de pós-processamento, análise de convergência e visualização gráfica em **Python**.

---

## Estrutura do Repositório

O repositório está organizado em Projetos de Programação Computacional (PPCs) e atividades práticas realizadas em sala de aula:

```plaintext
.
├── 📁 Atividades de sala                # Exercícios práticos e algoritmos fundamentais
│   ├── 📁 Distribuição de calor em placa # Condução bidimensional permanente
│   ├── 📁 EDO2_heun_euler               # Comparação de métodos de integração para EDOs de 2ª ordem
│   ├── 📁 Problema gravitacional de 2 corpos # Órbitas com Euler explícito vs. Leapfrog
│   ├── aurea.cpp                        # Otimização via Seção Áurea (1D)
│   ├── decomp_LU.cpp                    # Decomposição LU manual (10x10)
│   ├── interpolacao.cpp                 # Otimização por Interpolação Quadrática (1D)
│   ├── muller.cpp                       # Método de Müller para raízes complexas
│   ├── newton_raphson.cpp               # Método de Newton-Raphson para sistemas não-lineares
│   └── EDO_metodo_do_tiro.cpp           # Demonstração básica do Método do Tiro
│
├── 📁 PPC1 - Runge Kutta 4              # Sedimentação de partícula em regime Stokes/não-linear
├── 📁 PPC2 - Método de Bairstow         # Raízes de polinômios e geração de Fractais de Bairstow
├── 📁 PPC3 - Sistemas lineares          # Equação do calor 1D transiente via Thomas (Matriz Tridiagonal)
├── 📁 PPC4 - Otimização                 # Otimização 2D: Aclive Máximo vs. Gradientes Conjugados
├── 📁 PPC5 - Integradores numéricos     # Camada limite de Blasius via Método do Tiro + RK4
└── 📁 PPC6 - EDPs e Metodo de Liebmann  # Condução 2D em aleta: Gauss direto vs. Liebmann (Gauss-Seidel) + SOR
```

---

## Detalhamento dos Projetos de Programação (PPCs)

Cada diretório de PPC possui seu próprio `README.md` contendo a fundamentação teórica, formulações matemáticas e análises físicas detalhadas. Abaixo está um resumo de cada projeto:

| Projeto | Problema Físico / Matemático | Métodos Numéricos Utilizados | Destaques |
| :--- | :--- | :--- | :--- |
| **[PPC1](./PPC1%20-%20Runge%20Kutta%204)** | Sedimentação de uma partícula esférica em fluido (regimes de Stokes e inercial/Oseen). | **Runge-Kutta de 4ª Ordem (RK4)**. | Análise de erro global ($E \propto \Delta t^4$) e curvas de relaxação da velocidade. |
| **[PPC2](./PPC2%20-%20M%C3%A9todo%20de%20Bairstow)** | Busca de raízes reais e complexas para polinômios de ordem qualquer. | **Método de Bairstow** e **Newton-Raphson 2D**. | Geração de Fractais de Convergência de Bairstow no plano complexo. |
| **[PPC3](./PPC3%20-%20Sistemas%20lineares)** | Condução transiente de calor unidimensional em uma barra com geração de calor. | **Método das Diferenças Finitas (MDF)** Implícito + **Algoritmo de Thomas**. | Análise de estabilidade incondicional e animações temporais do perfil de temperatura. |
| **[PPC4](./PPC4%20-%20Otimiza%C3%A7%C3%A3o)** | Otimização multidimensional irrestrita (maximização de função quadrática 2D). | **Aclive Máximo** vs. **Gradientes Conjugados (Fletcher-Reeves)** + Busca Unidimensional por Interpolação Quadrática. | Comparação de caminhos de busca sobre curvas de nível e eficiência computacional. |
| **[PPC5](./PPC5%20-%20Integradores%20num%C3%A9ricos%20e%20EDOs)** | Escoamento sobre placa plana: Equação da Camada Limite de Blasius (PVC 3ª ordem). | **Método do Tiro (Shooting Method)** + Integração por **RK4**. | Determinação do fator de corte de Blasius e plotagem dos perfis de similaridade de velocidade. |
| **[PPC6](./PPC6%20-%20EDPs%20e%20Metodo%20de%20Liebmann)** | Distribuição bidimensional de temperatura em uma aleta retangular sob convecção. | **MDF 2D** + **Eliminação de Gauss** (Direto) vs. **Método de Liebmann / Gauss-Seidel** (Iterativo) com **SOR (Sobre-relaxação)**. | Estudo da influência do fator de sobre-relaxação ($\omega_{opt}$) e validação com solução analítica 1D. |

---

## Resumo das Atividades de Sala

O diretório **[Atividades de sala](./Atividades%20de%20sala)** reúne implementações eficientes e compactas de algoritmos clássicos estudados em aula:

1. **[aurea.cpp](./Atividades%20de%20sala/aurea.cpp)**: Otimização unidimensional para determinar a resistência elétrica ótima que maximiza a potência dissipada.
2. **[decomp_LU.cpp](./Atividades%20de%20sala/decomp_LU.cpp)**: Decomposição $LU$ e substituições diretas programadas sem bibliotecas externas para sistemas $10 \times 10$.
3. **[interpolacao.cpp](./Atividades%20de%20sala/interpolacao.cpp)**: Otimização unidimensional por interpolação parabólica (3 pontos).
4. **[muller.cpp](./Atividades%20de%20sala/muller.cpp)**: Método de Müller para cálculo de raízes complexas de polinômios através de parábolas aproximantes.
5. **[newton_raphson.cpp](./Atividades%20de%20sala/newton_raphson.cpp)**: Localização de raízes e resolução de sistemas não-lineares via Newton-Raphson.
6. **[Distribuição de calor em placa](./Atividades%20de%20sala/Distribui%C3%A7%C3%A3o%20de%20calor%20em%20placa)**: Método de diferenças finitas para condução 2D estacionária.
7. **[EDO2_heun_euler](./Atividades%20de%20sala/EDO2_heun_euler)**: Resolução de EDOs de 2ª ordem comparando Euler Simples e o método preditor-corretor de Heun.
8. **[Problema gravitacional de 2 corpos](./Atividades%20de%20sala/Problema%20gravitacional%20de%202%20corpos)**: Simulação orbital comparando a conservação de energia no método de Euler contra o método simplético Leapfrog.

---

## Pré-requisitos & Compilação

Para compilar e executar os códigos contidos neste repositório, você precisará de compiladores de C++ modernos e uma instalação funcional de Python para a geração de gráficos.

### C++

* Compilador compatível com **C++17** ou superior (recomenda-se `g++` 9.0+ ou `clang`).
* Makefile ou comandos diretos via terminal para compilação rápida.

### Python (Pós-processamento e Plotagem)

* Python 3.8+ instalado.
* Bibliotecas necessárias (instale via `pip`):

  ```bash
  pip install numpy matplotlib pandas seaborn polars science plots
  ```

> Nota: nem todos os pacotes são obrigatórios, os scrips podem ser facilmente alterados para usar somente `pandas` (`polars` não seria necessário) e para não precisar do scienceplots (comentar a linha `plt.style.use("science")`)

---

## Como Executar

Cada subpasta contém scripts orquestradores estruturados. De forma geral, o fluxo de execução segue as etapas:

### 1. Compilar o motor numérico em C++

Navegue até a pasta do projeto correspondente e execute:

```bash
# Exemplo para o PPC1
g++ -O3 -std=c++17 main.cpp -o main
```

> **Nota**: O uso da flag de otimização `-O3` é altamente recomendado para reduzir os tempos de processamento numérico, especialmente nas simulações mais densas como a geração de fractais no PPC2 ou malhas refinadas no PPC6.

### 2. Rodar o executável

```bash
./main
```

Isso gerará arquivos de dados (geralmente `.dat`, `.csv` ou `.txt`) contendo os resultados das simulações numéricas.

### 3. Executar o script de visualização em Python

```bash
python3 plot.py
# ou o respectivo arquivo de plotagem contido na pasta, ex:
python3 analysis.py
```

Esse comando lerá os dados salvos pelo resolvedor em C++ e criará os gráficos explicativos (salvos na pasta `figures/` ou no próprio diretório de trabalho).

---

## Sobre a Disciplina

* **Universidade**: Universidade de Brasília (UnB)
* **Faculdade**: Faculdade de Tecnologia (FT)
* **Departamento**: Departamento de Engenharia Mecânica (ENM)
* **Curso**: Engenharia Mecânica
* **Matéria**: Cálculo Numérico Aplicado
* **Professor**: Rafael Gabler Gontijo

---

## 📄 Licença

Este repositório está sob a licença MIT. Sinta-se livre para usar, estudar e modificar as implementações.
