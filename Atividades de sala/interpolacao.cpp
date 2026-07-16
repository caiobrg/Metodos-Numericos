// ============================================================================
// Método de Interpolação Quadrática para Otimização
// ----------------------------------------------------------------------------
// Encontra o valor de resistência ajustável (Ra) que maximiza a potência
// dissipada em um circuito elétrico. A cada iteração, ajusta-se uma parábola
// passando por 3 pontos e calcula-se o vértice como nova estimativa do ótimo.
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
  // --- Parâmetros do circuito elétrico ---
  double v = 80.0, r1 = 8.0, r2 = 12.0, r3 = 10.0;

  // --- Parâmetros de convergência ---
  double tol = 1.0e-4, erro = 1.0;

  // --- Três estimativas iniciais de Ra para a interpolação ---
  double ra_0 = 10.0, ra_1 = 39.0, ra_2 = 91.0, ra_3, ra_anterior;

  // Potências correspondentes a cada estimativa
  double p_0, p_1, p_2, p_3;

  // Proteção contra Ra = 0 (divisão por zero na função de potência)
  if (ra_0 == 0.0)
    ra_0 += 1.0e-4;

  // Lambda que calcula a abscissa do vértice da parábola interpoladora.
  // Dados 3 pontos (r0, p0), (r1, p1), (r2, p2), determina o Ra que
  // maximiza o polinômio quadrático ajustado via fórmula analítica do vértice.
  auto calculate_ra3 = [&](double r0, double r1, double r2) {
    double p0 = p(r0, v, r1, r2, r3);
    double p1 = p(r1, v, r1, r2, r3);
    double p2 = p(r2, v, r1, r2, r3);
    double num = p0 * (std::pow(r1, 2.0) - std::pow(r2, 2.0)) +
                 p1 * (std::pow(r2, 2.0) - std::pow(r0, 2.0)) +
                 p2 * (std::pow(r0, 2.0) - std::pow(r1, 2.0));
    double den =
        p0 * (r1 - r2) * 2.0 + p1 * (r2 - r0) * 2.0 + p2 * (r0 - r1) * 2.0;
    return num / den;
  };

  // --- Cálculo inicial: primeira estimativa ra_3 e potências nos 4 pontos ---
  ra_3 = calculate_ra3(ra_0, ra_1, ra_2);
  p_0 = p(ra_0, v, r1, r2, r3);
  p_1 = p(ra_1, v, r1, r2, r3);
  p_2 = p(ra_2, v, r1, r2, r3);
  p_3 = p(ra_3, v, r1, r2, r3);

  // --- Laço iterativo de otimização ---
  // A cada iteração, seleciona os 3 melhores pontos (lower, mid, upper)
  // que mantêm o máximo cercado, recalcula a parábola e atualiza ra_3.
  while (std::abs(erro) > tol) {
    ra_anterior = ra_3;
    double ra_lower, ra_mid, ra_upper;

    // Lógica de seleção dos 3 pontos para a próxima iteração:
    // Compara a posição e o valor de potência de ra_3 em relação a ra_1
    // para decidir qual dos 4 pontos descartar, mantendo o máximo confinado.
    if (ra_3 > ra_1) {
      if (p_3 > p_1) {
        // ra_3 está à direita de ra_1 e tem potência maior: descarta ra_0
        ra_lower = ra_1;
        ra_mid = ra_3;
        ra_upper = ra_2;
      } else {
        // ra_3 está à direita mas potência menor: descarta ra_2
        ra_lower = ra_0;
        ra_mid = ra_1;
        ra_upper = ra_3;
      }
    } else {
      if (p_3 > p_1) {
        // ra_3 está à esquerda de ra_1 e tem potência maior: descarta ra_2
        ra_lower = ra_0;
        ra_mid = ra_3;
        ra_upper = ra_1;
      } else {
        // ra_3 está à esquerda e potência menor: descarta ra_0
        ra_lower = ra_3;
        ra_mid = ra_1;
        ra_upper = ra_2;
      }
    }

    // Atualiza os 3 pontos e recalcula potências
    ra_0 = ra_lower;
    ra_1 = ra_mid;
    ra_2 = ra_upper;
    p_0 = p(ra_0, v, r1, r2, r3);
    p_1 = p(ra_1, v, r1, r2, r3);
    p_2 = p(ra_2, v, r1, r2, r3);

    // Calcula nova estimativa pelo vértice da parábola
    ra_3 = calculate_ra3(ra_0, ra_1, ra_2);
    p_3 = p(ra_3, v, r1, r2, r3);

    // Erro relativo entre estimativas consecutivas de Ra
    erro = (ra_3 - ra_anterior) / ra_3;
  }

  // --- Exibe resultado final ---
  std::cout << std::fixed << std::setprecision(8) << "Ra_3: " << ra_3
            << " | P_3: " << p_3 << std::endl;

  return 0;
}
