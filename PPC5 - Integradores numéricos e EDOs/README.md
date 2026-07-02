# PPC5/APC5 - Método do Tiro e Equação de Blasius

## Resumo

Este projeto implementa um solver numérico para a **equação de Blasius**, que descreve o perfil de velocidades da camada limite laminar de um escoamento incompressível sobre uma placa plana. O problema é um Problema de Valor de Contorno (PVC) não linear de terceira ordem, resolvido por meio da combinação de duas estratégias numéricas: o **Método do Tiro** (transformando o PVC em uma sequência de Problemas de Valor Inicial) e o **Método de Runge-Kutta de quarta ordem (RK4)** para a integração do sistema de EDOs de primeira ordem resultante. O projeto foi estruturado da seguinte maneira:

```plain text
└── 📁PPC5 - Metodo do Tiro
    ├── .gitignore
    ├── main.cpp
    ├── plot.py
    ├── README.md
    └── 📁results
        ├── iter1.csv
        ├── iter2.csv
        ├── ...
        └── iterN.csv
```

Todo o motor numérico (configuração, integrador RK4 e método do tiro) está implementado em um único arquivo orquestrador, [main.cpp](main.cpp), organizado em três classes independentes. A plotagem dos perfis de similaridade é feita pelo script em Python [plot.py](plot.py), que lê todos os arquivos gerados em `results/` e gera a figura `perfis_de_similaridade.png`.

Caso deseje somente as instruções de compilação e uso, favor se dirigir à seção [Instruções de Uso](#instruções-de-uso).

## 1 Introdução

O escoamento laminar bidimensional de um fluido incompressível sobre uma placa plana admite, por meio da variável de similaridade

$$\eta = y\sqrt{\dfrac{U_\infty}{\nu x}}$$

e da função de corrente $\psi = \sqrt{\nu U_\infty x}\, f(\eta)$, uma redução das equações da camada limite de Prandtl a uma única Equação Diferencial Ordinária não linear de terceira ordem, conhecida como **equação de Blasius**:

$$f''' + \frac{1}{2}ff'' = 0$$

sujeita às condições de contorno

$$f(0) = 0, \qquad f'(0) = 0, \qquad f'(\infty) = 1.$$

Como as condições de contorno estão distribuídas entre dois pontos distintos do domínio ($\eta = 0$ e $\eta \to \infty$), este é um PVC clássico. A estratégia adotada consiste em transformá-lo em um PVI através do **Método do Tiro**: o valor desconhecido $f''(0)$ é tratado como um parâmetro de ajuste $s$, integra-se o sistema por RK4 até um valor suficientemente grande de $\eta$ (aqui chamado $\eta_{max}$, em substituição a $\infty$) e corrige-se sucessivamente o valor de $s$ até que a condição $f'(\eta_{max}) \approx 1$ seja satisfeita dentro de uma tolerância especificada.

## 2 Formulação matemática

Reescrevendo a equação de Blasius como um sistema de três EDOs de primeira ordem, definem-se as variáveis auxiliares:

$$y_1 = f, \qquad y_2 = f', \qquad y_3 = f''$$

de modo que o sistema assume a forma

$$\frac{dy_1}{d\eta} = y_2, \qquad \frac{dy_2}{d\eta} = y_3, \qquad \frac{dy_3}{d\eta} = -\frac{1}{2}y_1 y_3$$

As condições de contorno originais são reescritas em termos dessas variáveis como $y_1(0) = 0$, $y_2(0) = 0$ e $y_2(\eta_{max}) = 1$. O parâmetro de tiro é justamente o valor inicial desconhecido $s = y_3(0) = f''(0)$: quanto maior o chute de $s$, mais rapidamente a solução tende à condição de contorno distante da parede, e o objetivo do método é encontrar o valor de $s$ que faz a solução chegar à condição de contorno $y_2(\eta_{max}) = 1$.

A função erro utilizada no Método do Tiro é $e(s) = y_2(\eta_{max}; s) - 1$, e a atualização dos sucessivos valores de $s$ é feita pelo **Método da Secante**, evitando a necessidade de calcular explicitamente a derivada de $e(s)$:

$$s_{i+1} = s_i - e(s_i)\frac{s_i - s_{i-1}}{e(s_i) - e(s_{i-1})}$$

## 3 Implementação

O código foi estruturado de forma orientada a objetos, em três classes com responsabilidades bem definidas:

* **`ConfiguracaoSimulacao`**: encapsula os parâmetros numéricos da simulação (chute inicial $s$, passo de integração $\Delta\eta$, $\eta_{max}$, tolerância e número máximo de iterações do tiro) e é responsável por capturá-los interativamente pelo terminal, com valores padrão caso o usuário não informe nada.
* **`IntegradorRK4`**: recebe um estado $(y_1, y_2, y_3)$ e avança um passo de integração $\Delta\eta$ utilizando o método de Runge-Kutta de quarta ordem, calculando os coeficientes $k_1$, $k_2$, $k_3$ e $k_4$ para as três variáveis do sistema simultaneamente.
* **`MetodoDoTiro`**: orquestra o processo iterativo completo. Realiza dois chutes iniciais para $s$, integra ambos até $\eta_{max}$, calcula os erros correspondentes e aplica o método da secante para gerar sucessivos valores até a convergência (ou até atingir o número máximo de iterações). A cada iteração, o perfil completo $(\eta, f, f', f'')$ é salvo em um arquivo `.csv` dentro de `results/`, permitindo visualizar a evolução da solução ao longo do processo de convergência.

Após a convergência, o programa realiza duas análises físicas adicionais a partir do perfil obtido:

1. **Coeficiente de atrito local**, isolando a constante $C_f\sqrt{Re_x} = 2f''(0)$;
2. **Espessura da camada limite** $\eta_{99}$, obtida por interpolação linear entre os dois pontos da malha que envolvem a condição $f'(\eta) = 0.99$, e comparada com a correlação clássica $\delta/x = 4.92/\sqrt{Re_x}$.

## 4 Resultados

Executando o programa com os parâmetros padrão ($s_0 = 0.2$, $\Delta\eta = 0.01$, $\eta_{max} = 10$, $tol = 10^{-5}$, $itMax = 100$), o Método do Tiro converge em apenas **6 iterações**, alternando entre os dois chutes iniciais e refinamentos sucessivos via secante:

| Iter | $f''(0) \rightarrow s$ | $f'(\eta_{max})$ | Erro |
| ------ | ----------- | ------------------- | ------ |
| 1 | 0.200000 | 0.713199 | 0.286801 |
| 2 | 0.400000 | 1.132134 | 0.132134 |
| 3 | 0.336919 | 1.009737 | 0.009737 |
| 4 | 0.331901 | 0.999686 | 0.000314 |
| 5 | 0.332058 | 1.000001 | 0.000001 |
| 6 | 0.332057 | 1.000000 | 0.000000 |

O valor convergido para $f''(0)$ foi de **0.332057337234546**, em excelente concordância com o valor clássico da literatura $f''(0) \approx 0.332057$, com erro final $|f'(\eta_{max}) - 1| \approx 6.0\times10^{-11}$.

A partir desse valor, obteve-se a constante do coeficiente de atrito local $C_f\sqrt{Re_x} = 0.664114674469091$ e a posição adimensional $\eta_{99} = 4.909989482810054$, correspondente ao coeficiente de espessura $C_\delta = \eta_{99} \approx 4.910$. Esse resultado está muito próximo da correlação clássica $C_\delta = 4.92$, com uma diferença relativa inferior a $0.2\%$, o que valida tanto a implementação do integrador RK4 quanto a estratégia do Método do Tiro.

As principais fontes de erro numérico residual estão associadas:

(i) ao passo de integração $\Delta\eta$, cujo refinamento reduz o erro de truncamento local do RK4 (erro global da ordem $\Delta\eta^4$);
(ii) à escolha de $\eta_{max}$, que precisa ser suficientemente grande para que a condição $f'(\eta) \to 1$ seja de fato assintótica e não apenas aproximada dentro do domínio truncado;
(iii) à tolerância adotada no Método do Tiro, que define o critério de parada da secante.

O script `plot.py` gera a figura `perfis_de_similaridade.png`, com os perfis de $f(\eta)$, $f'(\eta)$ e $f''(\eta)$ obtidos em cada iteração do método do tiro, permitindo visualizar tanto a solução final convergida quanto a evolução dos chutes intermediários.

## Instruções de Uso

### Pré-requisitos

* **Compilador C++**: GCC (`g++`) ou Clang com suporte a C++17 ou superior (devido ao uso da biblioteca `<filesystem>`).
* **Python 3.x**: recomenda-se a instalação via `pip` das bibliotecas `polars`, `matplotlib` e `scienceplots`:

```bash
pip install polars matplotlib scienceplots
```

> Observação: o estilo `science` do `scienceplots` depende de uma instalação de LaTeX no sistema. Caso o LaTeX não esteja disponível, remova ou comente a linha `plt.style.use("science")` em `plot.py` para gerar a figura com o estilo padrão do Matplotlib.

### Passo 1: Compilação do Orquestrador (C++)

Compile o código ativando otimizações:

```bash
g++ -O3 -std=c++17 main.cpp -o main
```

### Passo 2: Execução

Ao executar o programa, ele solicitará interativamente os parâmetros numéricos da simulação. Basta pressionar Enter para aceitar os valores padrão indicados entre colchetes:

```bash
./main
```

* `s` (chute inicial de $f''(0)$): default `0.2`;
* `\delta\eta` (passo de integração): default `0.01`;
* `\eta_{max}` (ponto máximo da integração): default `10.0`;
* `tol` (tolerância de convergência): default `1e-5`;
* `itMax` (máximo de iterações do Método do Tiro): default `100`.

O programa criará automaticamente o diretório `results/`, salvará o perfil de cada iteração do tiro em arquivos `iterN.csv` (colunas `eta,y1,y2,y3`) e imprimirá na tela:

* a evolução do chute, do valor de $f'(\eta_{max})$ e do erro a cada iteração;
* o valor convergido de $f''(0)$;
* o número de iterações necessárias para a convergência;
* o erro final obtido;
* o valor final de $f'(\eta_{max})$;
* a constante do coeficiente de atrito $C_f\sqrt{Re_x}$;
* a posição adimensional $\eta_{99}$ e o coeficiente de espessura $C_\delta$.

É possível alterar este comportamento alterando variáveis no orquestrador ([main.cpp](main.cpp)), sendo possível configurar o programa para salvar e escrever apenas os resultados finais.

### Passo 3: Plotagem dos Perfis de Similaridade

Com o diretório `results/` já populado, execute o script de plotagem no mesmo diretório:

```bash
python plot.py
```

A figura `perfis_de_similaridade.png` será salva no diretório raiz, contendo os perfis de $f(\eta)$, $f'(\eta)$ e $f''(\eta)$ para cada iteração do Método do Tiro.

### Passo 4: Verificação de Saídas

Os resultados poderão ser conferidos nos seguintes arquivos:

1. `results/iterN.csv`: perfis $(\eta, f, f', f'')$ de cada iteração do Método do Tiro, sendo o último arquivo gerado o correspondente à solução convergida;
2. `perfis_de_similaridade.png`: gráfico comparativo dos perfis de similaridade de Blasius ao longo das iterações do tiro.
