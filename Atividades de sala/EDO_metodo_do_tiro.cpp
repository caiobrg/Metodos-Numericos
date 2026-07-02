#include <iostream>
#include <ostream>
#include <vector>

int main() {
  // Parâmetros físicos
  double L = 1;       // m, comprimento da barra
  int N = 100;        // Número de pontos da malha
  double q_dot = 1e6; // W/m^3, Geração interna de calor
  double k = 60;      // W/(m*K), condutividade do material

  // Condições de contorno
  double T_0 = -1000;  // Temperatura no ponto 0
  double T_L = 100000; // Temperatura no ponto L

  double dx = L / (N - 1); // m, Espaçamento da malha

  // Construção do vetor de malha
  std::vector<double> malha;
  for (double x = 0; x <= L; x += dx) {
    malha.push_back(x);
  }

  // Variáveis de convergência
  double tol = 1e-6;
  double err1 = 1, err2 = 1;

  // Chutes iniciais e taxa de variação no espaço
  double df_dx = -q_dot / k;
  double f_1_ini = -100;
  double f_2_ini = 20;

  //   Variáveis auxiliares
  // T_i: temperaturas para os chutes 1 e 2 para integração pelo método de Euler
  // f_3: variável para atualização dos valores iniciais pelo método da secante
  std::vector<double> T_1(N, T_0);
  std::vector<double> T_2(N, T_0);
  double f_3;

  // Run inicial para obtenção de 2 chutes

  // f_i: taxas de variação da temperatura para os chutes 1 e 2
  std::vector<double> f_1(N, f_1_ini);
  std::vector<double> f_2(N, f_2_ini);

  for (int i = 1; i < malha.size(); i++) {
    T_1[i] = (T_1[i - 1] + f_1[i - 1] * dx);
    T_2[i] = (T_2[i - 1] + f_2[i - 1] * dx);
    f_1[i] = (f_1[i - 1] + df_dx * dx);
    f_2[i] = (f_2[i - 1] + df_dx * dx);
  }

  err1 = T_1[N - 1] - T_L;

  double iter = 0;

  while (std::abs(err1) >= tol or std::abs(err2) >= tol) {

    err2 = T_2[N - 1] - T_L;
    f_3 = f_2[0] - err2 * (f_2[0] - f_1[0]) / (err2 - err1);
    f_1[0] = f_2[0];
    err1 = err2;
    f_2[0] = f_3;

    for (int i = 1; i < malha.size(); i++) {
      T_2[i] = (T_2[i - 1] + f_2[i - 1] * dx);
      f_2[i] = (f_2[i - 1] + df_dx * dx);
    }
    iter++;
  }
  std::cout << err1 << " " << err2 << std::endl;
  std::cout << T_2[N - 1] << std::endl;
  std::cout << iter << std::endl;
  return 0;
}
