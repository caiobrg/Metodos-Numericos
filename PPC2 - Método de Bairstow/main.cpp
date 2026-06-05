#include "bairstow_solver.hpp"
#include <fstream> // Para salvar arquivos no computador
#include <iomanip> // Usado para deixar a impressão dos números mais bonita
#include <iostream>
#include <vector>

/**
 * @brief Realiza a validação do solver com um polinômio conhecido.
 *
 * @param polinomio Vetor de coeficientes do polinômio.
 * @param tolerance Tolerância para convergência.
 * @param maxIterations Limite de iterações.
 */
void validacao(std::vector<double> polinomio, double tolerance = 1e-6,
               int maxIterations = 1000) {
  std::cout << "========================================\n";
  std::cout << "       1. Validação do solver           \n";
  std::cout << "========================================\n\n";

  // Chutes iniciais para o fator quadrático x^2 - rx - s
  double rGuess = 2.5;
  double sGuess = 2.5;

  std::cout << "Palpites iniciais: r = " << rGuess << ", s = " << sGuess
            << "\n";
  std::cout << "Tolerancia: " << tolerance << "\n\n";

  // Chamada do algoritmo de Bairstow
  bairstowResult resultado =
      bairstowSolve(polinomio, rGuess, sGuess, tolerance, maxIterations);

  // Configuração da formatação de saída
  std::cout << std::fixed << std::setprecision(6);

  for (size_t i = 0; i < resultado.roots.size(); ++i) {
    std::cout << "Raiz " << i + 1 << ": ";
    std::cout << resultado.roots[i].real();

    // Exibição da parte imaginária com sinal adequado
    if (resultado.roots[i].imag() >= 0) {
      std::cout << " + " << resultado.roots[i].imag() << "i\n";
    } else {
      std::cout << " - " << std::abs(resultado.roots[i].imag()) << "i\n";
    }
  }
  std::cout << "----------------------------------------\n";

  std::cout << "Total de iterações: " << resultado.total_iterations << "\n";
  std::cout << "Processo convergiu: "
            << ((resultado.success == 1) ? "Sim" : "Não") << "\n";
}

/**
 * @brief Analisa a convergência do método para diferentes chutes iniciais.
 */
void analiseConvergencia(std::vector<double> polinomio, double tolerance = 1e-6,
                         int maxIterations = 1000) {
  std::cout << "========================================\n";
  std::cout << "     2. Análise de Convergência         \n";
  std::cout << "========================================\n\n";

  std::vector<double> rGuesses = {0.1, 0.25, 0.5, 1, 2, 5, 10};
  std::vector<double> sGuesses = {0.1, 0.25, 0.5, 1, 2, 5, 10};

  for (int i = 0; i < rGuesses.size(); i++) {
    bairstowResult resultado = bairstowSolve(
        polinomio, rGuesses[i], sGuesses[i], tolerance, maxIterations);

    std::cout << "Par (r, s) do chute (" << rGuesses[i] << ", " << sGuesses[i]
              << ")\n";
    std::cout << "Total de iterações: " << resultado.total_iterations << "\n";
    std::cout << "Processo convergiu: "
              << ((resultado.success == 1) ? "Sim" : "Não") << "\n";
  }
}

/**
 * @brief Aplicação do método a um caso específico de engenharia, um sistema
 * massa-mola-amortecedor.
 */
void analiseDeCaso(std::vector<double> polinomio, double tolerance = 1e-6,
                   int maxIterations = 1000) {
  std::cout << "==========================================\n";
  std::cout << "3. Analise de Caso - Massa Mola Amortecido\n";
  std::cout << "==========================================\n\n";

  double rGuess = 2.5;
  double sGuess = 2.5;

  std::cout << "Palpites iniciais: r = " << rGuess << ", s = " << sGuess
            << "\n";
  std::cout << "Tolerancia: " << tolerance << "\n\n";

  bairstowResult resultado =
      bairstowSolve(polinomio, rGuess, sGuess, tolerance, maxIterations);

  std::cout << std::fixed << std::setprecision(6);

  for (size_t i = 0; i < resultado.roots.size(); ++i) {
    std::cout << "Raiz " << i + 1 << ": ";
    std::cout << resultado.roots[i].real();

    if (resultado.roots[i].imag() >= 0) {
      std::cout << " + " << resultado.roots[i].imag() << "i\n";
    } else {
      std::cout << " - " << std::abs(resultado.roots[i].imag()) << "i\n";
    }
  }
  std::cout << "----------------------------------------\n";

  std::cout << "Total de iterações: " << resultado.total_iterations << "\n";
  std::cout << "Processo convergiu: "
            << ((resultado.success == 1) ? "Sim" : "Não") << "\n";
}

struct Range {
  double min;
  double max;
};

/**
 * @brief Gera os dados para a visualização do Fractal de Bairstow.
 *
 * Varre um plano de chutes iniciais (r, s) e registra a convergência.
 */
void fractalDeBairstow(std::vector<double> polinomio,
                       const char *filename = "fractal.csv",
                       Range rRange = {-3.0, 3.0}, Range sRange = {-3.0, 3.0},
                       int maxIterations = 100, int res = 5000,
                       double tolerance = 1e-6) {
  std::cout << "========================================\n";
  std::cout << "         4. Fractal de Bairstow         \n";
  std::cout << "========================================\n\n";

  std::cout << "Resultados serão salvos no arquivo: " << filename << "\n";

  std::ofstream arquivo_csv(filename);

  if (!arquivo_csv.is_open()) {
    std::cerr << "ERRO: Nao foi possivel criar o arquivo CSV" << "\n";
    return;
  }

  // Cabeçalho do CSV: r e s iniciais, iterações, se convergiu, e os valores de
  // r e s convergidos
  arquivo_csv << "r,s,iteracoes,sucesso,r_conv,s_conv,exact_root\n";

  double rStep = (rRange.max - rRange.min) / res;
  double sStep = (sRange.max - sRange.min) / res;

  std::cout << "Iniciando a varredura (Resolução: " << res << "x" << res
            << ")\n";

  // Paralelização com OpenMP para acelerar a varredura do plano
#pragma omp parallel for collapse(2)
  for (int i = 0; i <= res; i++) {
    for (int j = 0; j <= res; j++) {
      double r = rRange.min + (i * rStep);
      double s = sRange.min + (j * sStep);

      bairstowResult resultado =
          bairstowSolve(polinomio, r, s, tolerance, maxIterations);

      // Seção crítica para escrita segura no arquivo (evita race conditions)
#pragma omp critical
      {
        arquivo_csv << r << "," << s << "," << resultado.total_iterations << ","
                    << resultado.success << "," << resultado.r_converged << ","
                    << resultado.s_converged << "," << resultado.exact_root
                    << "\n";
      }
    }
  }
}

int main() {
  // Criação da pasta de outputs se não existir
  system("mkdir -p outputs");

  std::cout << "========================================\n";
  std::cout << "   Orquestrador - Metodo de Bairstow    \n";
  std::cout << "========================================\n\n";

  std::cout << "P(x) = (x-1)(x+2)(x-3)(x+4)(x-5)(x+6)(x-7) = x^7 - 4 "
               "x^6 - 62 x^5 + 200 x^4 + 1009 x^3 - 2356 x^2 - 3828 x + 5040\n";

  // Definindo a lista de coeficiêntes do polinomio
  std::vector<double> polinomio = {5040, -3828, -2356, 1009, 200, -62, -4, 1};

  double tolerance = 1e-6;
  int maxIterations = 1000;

  validacao(polinomio, tolerance, maxIterations);
  analiseConvergencia(polinomio, tolerance, maxIterations);

  // Polinômio característico (Grau 6)
  std::cout << "\nP(x) = 6 x^6 + 25/2 x^5 + 253/2 x^4 + 601/8 x^3 + 3201/8 "
               "x^2 + 45 x + 162\n";
  std::vector<double> polinomioCaract = {
      162, 45, 3201.0 / 8.0, 601.0 / 8.0, 253.0 / 2.0, 12.5, 6};

  analiseDeCaso(polinomioCaract);

  maxIterations = 100;
  int resolucao = 3000;
  fractalDeBairstow(polinomioCaract, "outputs/fractal_pol_carac.csv", {-10, 10},
                    {-10, 10}, maxIterations, resolucao);
  fractalDeBairstow(polinomio, "outputs/fractal_pol_7.csv", {-10, 10},
                    {-10, 10}, maxIterations, resolucao);
  fractalDeBairstow({53, 4, 23, -1067, -12, 8, -2},
                    "outputs/fractal_pol_comp.csv", {-5, 5}, {-5, 5},
                    maxIterations, resolucao);

  return 0;
}
