#pragma once
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <variant>
#include <vector>

inline void solver_metodo_thomas(
    std::variant<double, std::vector<double>> temp_ini, double temp_inf,
    double L, int n_pontos, double t_total, double dt, double beta,
    double geracao, double t_ini = 0, double coef_conv1 = 15,
    double coef_conv2 = 15, double condutiv = 400, double calor_espec = 470,
    double dens = 7800, bool verbose = false, bool save_results = false,
    const char *filename = "results.bin") {

  std::ofstream arquivo_bin(filename, std::ios::binary);

  double dx = L / (n_pontos - 1);
  double k = condutiv;
  double cp = calor_espec;
  double rho = dens;
  double alpha = k / (cp * rho);
  double h_1 = coef_conv1;
  double h_2 = coef_conv2;
  double q_dot = geracao;

  double Fo = (alpha * dt) / (dx * dx);
  double Bi_1 = (h_1 * dx) / k;
  double Bi_2 = (h_2 * dx) / k;
  double A = (q_dot * dt) / (rho * cp);

  double um_men_beta = (1.0 - beta);

  std::vector<double> T(n_pontos);

  if (verbose) {
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Parâmetros:\nDifusividade = " << alpha << "\nFo = " << Fo
              << "\nBi_1 = " << Bi_1 << "\nBi_2 = " << Bi_2 << "\nA = " << A
              << std::endl;
  }

  // Verifica se o t_ini é um vetor ou um double
  if (std::holds_alternative<double>(temp_ini)) {
    // Se for um double preenche o vetor de temperaturas com ele
    double valor_constante = std::get<double>(temp_ini);
    T.assign(n_pontos, valor_constante);
  } else {
    // Se for um vetor apenas copia o vetor passado para o vetor T
    T = std::get<std::vector<double>>(temp_ini);

    // Se o vetor de temperaturas fornecido não possuir o mesmo número de
    // pontos n_pontos levanta um erro.
    if (T.size() != n_pontos) {
      throw "Erro: tamanho do vetor inicial não bate com o tamanho da malha";
    }
  }

  if (verbose)
    std::cout << "Iniciando solução do sistema linear" << std::endl;
  // Cria os vetores da matriz preenchidos com os valores que aparecem mais
  // vezes
  std::vector<double> e(n_pontos, -Fo * um_men_beta);
  std::vector<double> f(n_pontos, 1 + 2 * Fo * um_men_beta);
  std::vector<double> g(n_pontos, -Fo * um_men_beta);

  // Definindo os valores das extremidades dos vetores
  e[n_pontos - 1] = g[0] = -2 * Fo * um_men_beta;
  f[0] = 1 + (2 * Fo) * um_men_beta + (2 * Bi_1 * Fo) * um_men_beta;
  f[n_pontos - 1] = 1 + (2 * Fo) * um_men_beta + (2 * Bi_2 * Fo) * um_men_beta;

  // Vetores upper e lower a serem calculados no próximo loop
  std::vector<double> u(n_pontos, 0);
  std::vector<double> l(n_pontos, 0);

  // Valor inicial que não é calculado no loop
  u[0] = f[0];

  for (size_t k = 1; k < n_pontos; k++) {
    l[k] = e[k] / u[k - 1];
    u[k] = f[k] - (l[k] * g[k - 1]);
  }

  std::vector<double> b(n_pontos, 0);

  double tempo = t_ini;
  int passo_atual = 0;

  while (tempo < t_total) {
    // Definindo os valores das extremidades
    b[0] = T[0] + 2 * Fo * beta * (T[1] - T[0] - Bi_1 * T[0]) +
           2 * Fo * Bi_1 * temp_inf + A;
    b[n_pontos - 1] =
        T[n_pontos - 1] +
        2 * Fo * beta *
            (T[n_pontos - 2] - T[n_pontos - 1] - Bi_2 * T[n_pontos - 1]) +
        2 * Fo * Bi_2 * temp_inf + A;

    // Definindo os valores "internos"
    for (size_t m = 1; m < n_pontos - 1; m++) {
      b[m] = A + T[m] - Fo * beta * (2 * T[m] - T[m - 1] - T[m + 1]);
    }

    std::vector<double> d(n_pontos, 0);
    // Valor não calculado no loop
    d[0] = b[0];

    // Substituição progressiva
    for (size_t k = 1; k < n_pontos; k++) {
      d[k] = b[k] - l[k] * d[k - 1];
    }

    std::vector<double> x(n_pontos, 0);
    // Valor não calculado no loop
    x[n_pontos - 1] = d[n_pontos - 1] / u[n_pontos - 1];

    // Substituição regressiva
    for (int k = n_pontos - 2; k >= 0; k--) {
      x[k] = (d[k] - (g[k] * x[k + 1])) / u[k];
    }

    if (save_results)
      arquivo_bin.write(reinterpret_cast<const char *>(T.data()),
                        n_pontos * sizeof(double));

    // O vetor de temperaturas x vira o novo vetor de "entrada" de temperaturas
    T = x;

    // Atualiza o passo de tempo
    tempo += dt;
    passo_atual++;
  }
  if (verbose)
    std::cout << "Programa finalizado, dados gravados em: " << filename
              << std::endl;
}