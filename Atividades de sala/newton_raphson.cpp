// ============================================================================
// Método de Newton-Raphson para Localização de Raízes
// ----------------------------------------------------------------------------
// Encontra raízes de equações algébricas não lineares usando o método de
// aproximações sucessivas de Newton-Raphson clássico.
// Inclui variante modificada para tratamento de raízes múltiplas.
// Função alvo: f(x) = x³ - 6x² + 9x - 4
// ============================================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

// Função polinomial cuja raiz se deseja encontrar: f(x) = x³ - 6x² + 9x - 4
double f(double x) {
  return std::pow(x, 3.0) - 6.0 * std::pow(x, 2.0) + 9.0 * x - 4.0;
}

// Primeira derivada: f'(x) = 3x² - 12x + 9
double f_derivada1(double x) { return 3.0 * std::pow(x, 2.0) - 12.0 * x + 9.0; }

// Segunda derivada: f''(x) = 6x - 12
double f_derivada2(double x) { return 6.0 * x - 12.0; }

// Newton-Raphson clássico: x_{n+1} = x_n - f(x_n) / f'(x_n)
double f_NR(double x) { return x - (f(x) / f_derivada1(x)); }

// Newton-Raphson modificado para raízes múltiplas:
// x_{n+1} = x_n - [f(x)*f'(x)] / [f'(x)² - f(x)*f''(x)]
// Converge quadraticamente mesmo para raízes com multiplicidade > 1
double f_NRmod(double x) {
  return x - ((f(x) * f_derivada1(x)) /
              (std::pow(f_derivada1(x), 2.0) - (f(x) * f_derivada2(x))));
}

int main() {
  // --- Configuração inicial ---
  double x_atual = 100;       // Chute inicial
  double tol = 1.0e-6;        // Tolerância de convergência
  int max_iter = 100, iter = 0;
  double erro = 1.0;          // Erro absoluto |x_{n+1} - x_n|

  // Cabeçalho da tabela de saída
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << "                       METODO DE NEWTON-RAPHSON                "
               "                 "
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << " Chute inicial: " << std::fixed << std::setprecision(4)
            << x_atual << std::endl;
  std::cout << " Tolerancia:    " << std::scientific << tol << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << std::setw(6) << " Iter " << " | " << std::setw(15)
            << "    x_atual     " << " | " << std::setw(17)
            << "      f(x)       "
            << " | " << std::setw(17) << "      Erro       " << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  // --- Laço iterativo de Newton-Raphson ---
  while (erro > tol && iter < max_iter) {
    iter++;
    double fx = f(x_atual);         // Valor da função no ponto atual
    double df = f_derivada1(x_atual); // Derivada no ponto atual

    // Calcula a próxima estimativa usando NR clássico (ou NR modificado)
    double x_novo = f_NR(x_atual); // or f_NRmod(x_atual)
    erro = std::abs(x_novo - x_atual);

    // Exibe resultados da iteração
    std::cout << std::setw(5) << iter << " | " << std::fixed
              << std::setprecision(8) << std::setw(15) << x_atual << " | "
              << std::scientific << std::setprecision(8) << std::setw(17) << fx
              << " | " << std::setw(17) << erro << std::endl;

    // Proteção contra derivada nula (divisão por zero)
    if (std::abs(df) < 1e-12) {
      std::cout << "Aviso: Derivada próxima de zero. Parando." << std::endl;
      break;
    }

    x_atual = x_novo; // Atualiza a estimativa
  }

  // --- Resultado final ---
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << " Raiz encontrada:      " << std::fixed << std::setprecision(8)
            << x_atual << std::endl;
  std::cout << " Iteracoes realizadas: " << iter << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  return 0;
}
