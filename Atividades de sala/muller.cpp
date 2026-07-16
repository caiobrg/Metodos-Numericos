// ============================================================================
// Método de Muller para Localização de Raízes
// ----------------------------------------------------------------------------
// Generalização do método da secante que utiliza interpolação quadrática
// (parábola) através de 3 pontos para projetar o cruzamento com o eixo x.
// Aplicado a um polinômio de grau 3: f(x) = x³ - 13x - 12.
// Resultado salvo no terminal e no arquivo "raizes.dat".
// ============================================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// Função polinomial cuja raiz se deseja encontrar: f(x) = x³ - 13x - 12
double f(double x) { return std::pow(x, 3.0) - 13.0 * x - 12.0; }

int main() {
  // --- Chutes iniciais igualmente espaçados ---
  double x0 = 6.5, passo = 0.5;
  double x1 = x0 - passo;   // x1 = 6.0
  double x2 = x1 - passo;   // x2 = 5.5
  double x3 = 0.0;           // Nova estimativa (calculada a cada iteração)

  // --- Parâmetros de convergência ---
  double tol = 1.0e-6;       // Tolerância absoluta
  double erro = 1.0;         // Erro atual |x3 - x2|
  int iter = 0, max_iter = 100;

  // Cabeçalho da tabela de saída
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << "                              METODO DE MULLER                 "
               "                 "
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::fixed << std::setprecision(4)
            << " Chutes iniciais (x0, x1, x2): " << x0 << ", " << x1 << ", "
            << x2 << std::endl;
  std::cout << " Tolerancia:                   " << std::scientific << tol
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << std::setw(6) << " Iter " << " | " << std::setw(15)
            << "      x3       "
            << " | " << std::setw(17) << "      f(x3)      " << " | "
            << std::setw(17) << "      Erro       " << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  // Abre arquivo para salvar a raiz encontrada
  std::ofstream outfile("raizes.dat");

  // --- Laço iterativo do Método de Muller ---
  while (erro > tol && iter < max_iter) {
    iter++;

    // Cálculo dos passos entre os pontos consecutivos
    double h0 = x1 - x0; // Passo entre x0 e x1
    double h1 = x2 - x1; // Passo entre x1 e x2

    // Diferenças divididas de primeira ordem
    double delta0 = (f(x1) - f(x0)) / h0;
    double delta1 = (f(x2) - f(x1)) / h1;

    // Coeficientes do polinômio quadrático de ajuste local: a*y² + b*y + c = 0
    double a = (delta1 - delta0) / (h1 + h0); // Coef. quadrático (diferença dividida de 2ª ordem)
    double b = a * h1 + delta1;                // Coef. linear
    double c = f(x2);                           // Coef. independente (valor da função em x2)

    // Discriminante da equação quadrática
    double discriminante = std::pow(b, 2.0) - 4.0 * a * c;

    // Verificação de raízes complexas (discriminante negativo)
    if (discriminante < 0.0) {
      std::cout << "Aviso: Discriminante negativo. Raízes complexas detectadas."
                << std::endl;
      break;
    }

    // Escolha do sinal do denominador para maximizar |den|,
    // garantindo o menor passo possível e maior estabilidade numérica
    double den;
    if (b >= 0.0) {
      den = b + std::sqrt(discriminante);
    } else {
      den = b - std::sqrt(discriminante);
    }

    // Calcula a nova estimativa da raiz
    x3 = x2 - (2.0 * c) / den;
    erro = std::abs(x3 - x2);

    // Exibe resultados da iteração
    std::cout << std::setw(5) << iter << " | " << std::fixed
              << std::setprecision(8) << std::setw(15) << x3 << " | "
              << std::scientific << std::setprecision(8) << std::setw(17)
              << f(x3) << " | " << std::setw(17) << erro << std::endl;

    // Desloca os 3 pontos para a próxima iteração
    x0 = x1;
    x1 = x2;
    x2 = x3;
  }

  // --- Resultado final ---
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << " Raiz final encontrada: " << std::fixed << std::setprecision(8)
            << x3 << std::endl;
  std::cout << " Iteracoes realizadas:  " << iter << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  // Salva a raiz no arquivo de saída
  if (outfile.is_open()) {
    outfile << "Raiz encontrada: " << std::fixed << std::setprecision(8) << x3
            << std::endl;
    outfile.close();
  }

  return 0;
}
