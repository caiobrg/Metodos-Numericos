// ============================================================================
// Método da Seção Áurea para Otimização Unidimensional
// ----------------------------------------------------------------------------
// Encontra o valor de resistência ajustável (Ra) que maximiza a potência
// dissipada em um circuito elétrico. Utiliza a razão áurea para reduzir
// iterativamente o intervalo de busca até atingir a tolerância desejada.
// ============================================================================

#include <cmath>
#include <iomanip>
#include <iostream>

// Função de potência dissipada no circuito
// P(Ra) = (V * R3 * Ra / denominador)^2 / Ra
// Parâmetros do circuito: v = tensão, r1/r2/r3 = resistências fixas, ra = resistência ajustável
double p(double ra, double v, double r1, double r2, double r3) {
  double num = v * r3 * ra;
  double den = r1 * (ra + r2 + r3) + r3 * ra + r3 * r2;
  return std::pow(num / den, 2.0) / ra;
}

int main() {
  double ra_opt, p_max, ra1, ra2;

  // --- Parâmetros do circuito elétrico ---
  double r1 = 8.0, r2 = 12.0, r3 = 10.0;
  double v = 80.0; // Tensão da fonte [V]

  // Razão áurea: phi = (sqrt(5) - 1) / 2 ≈ 0.618033
  double razao = (std::sqrt(5.0) - 1.0) / 2.0;

  // --- Intervalo inicial de busca para Ra ---
  double ra_upper = 100.0; // Limite superior [Ω]
  double ra_lower = 0.0;   // Limite inferior [Ω]
  double length_1;          // Comprimento do subintervalo áureo

  // --- Parâmetros de convergência ---
  double tol = 1.0e-4; // Tolerância de parada
  double erro = 1.0;    // Erro atual (diferença entre ra1 e ra2)
  int iter = 0;         // Contador de iterações

  // Cabeçalho da tabela de saída
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << "                      OTIMIZACAO: METODO DA RAZAO AUREA        "
               "                 "
            << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;
  std::cout << std::setw(6) << " Iter " << " | " << std::setw(14)
            << "   ra_lower    " << " | " << std::setw(14) << "   ra_upper    "
            << " | " << std::setw(14) << "  ra_otimizado " << " | "
            << std::setw(14) << "      Erro      " << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-----------------"
            << std::endl;

  // --- Laço iterativo do método da seção áurea ---
  // A cada iteração, dois pontos internos (ra1, ra2) particionam o intervalo
  // usando a razão áurea. O subintervalo com menor potência é descartado,
  // reduzindo o intervalo por um fator de phi a cada passo.
  while (std::abs(erro) > tol) {
    iter++;

    // Calcula o comprimento do subintervalo pela razão áurea
    length_1 = (ra_upper - ra_lower) * razao;

    // Pontos internos de amostragem
    ra1 = ra_lower + length_1; // Ponto interno superior
    ra2 = ra_upper - length_1; // Ponto interno inferior

    // Compara potências nos pontos internos e descarta o subintervalo inferior
    if (p(ra2, v, r1, r2, r3) > p(ra1, v, r1, r2, r3)) {
      ra_upper = ra1; // Descarta a porção direita (ra1 vira novo limite superior)
    } else {
      ra_lower = ra2; // Descarta a porção esquerda (ra2 vira novo limite inferior)
    }

    // Atualiza o erro (distância entre pontos internos) e a melhor estimativa
    erro = ra1 - ra2;
    ra_opt = (ra_upper + ra_lower) / 2.0;

    // Exibe os dados da iteração atual
    std::cout << std::setw(6) << iter << " | " << std::fixed
              << std::setprecision(8) << std::setw(14) << ra_lower << " | "
              << std::setw(14) << ra_upper << " | " << std::setw(14) << ra_opt
              << " | " << std::scientific << std::setprecision(6)
              << std::setw(14) << erro << std::endl;
  }

  // --- Resultado final ---
  ra_opt = (ra_upper + ra_lower) / 2.0;
  p_max = p(ra_opt, v, r1, r2, r3);

  std::cout << "---------------------------------------------------------------"
               "--------------------"
            << std::endl;
  std::cout << "Ra otimizado: " << std::fixed << std::setprecision(8) << ra_opt
            << std::endl;
  std::cout << "Potência máxima: " << p_max << std::endl;
  std::cout << "---------------------------------------------------------------"
               "--------------------"
            << std::endl;

  return 0;
}
