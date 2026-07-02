#include <cmath>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

double dv_dt(double v, double x, double t, double m, double k, double c,
             double f_0, double w) {
  return (-k * x - c * v + f_0 * std::sin(w * t)) / m;
}

void integracao_euler(double t_total, double delta_t, double m, double c,
                      double k, double f_0, double w, double x_ini = 0.0,
                      double v_ini = 0.0, bool save_to_file = true,
                      const char *filename = "integração_euler.dat") {

  int N = std::round(t_total / delta_t + 0.5); // número de passos de simulação

  // Vetores temporais de x e v, inicializados completamente com os valores
  // iniciais
  std::vector<double> x(N, x_ini);
  std::vector<double> v(N, v_ini);

  for (int i = 1; i < N; i++) {
    double t = delta_t * i;
    v[i] = v[i - 1] +
           delta_t * dv_dt(v[i - 1], x[i - 1], (t - delta_t), m, k, c, f_0, w);
    x[i] = x[i - 1] + delta_t * v[i];
  }

  // Arquivo de saída
  if (save_to_file) {
    std::ofstream outfile(filename);

    if (outfile.is_open()) {
      outfile << "t,x,v" << std::endl;
      for (int i = 0; i < x.size(); i++) {
        outfile << i * delta_t << "," << x[i] << "," << v[i] << std::endl;
      }
      outfile.close();
      std::cout << "Simulação pelo métodod de Euler concluída com sucesso"
                << std::endl;
    } else {
      std::cout << "Não foi possível abrir o arquivo de saída" << std::endl;
    }
  } else {
    std::cout << "t,x,v" << std::endl;
    for (int i = 0; i < x.size(); i++) {
      std::cout << i * delta_t << "," << x[i] << "," << v[i] << std::endl;
    }
  }
}

void integracao_heun(double t_total, double delta_t, double m, double c,
                     double k, double f_0, double w, double x_ini = 0.0,
                     double v_ini = 0.0, bool save_to_file = true,
                     const char *filename = "integração_heun.dat") {
  int N = t_total / delta_t;
  std::vector<double> x(N, x_ini);
  std::vector<double> v(N, v_ini);
  std::vector<double> x_0(N, x_ini);
  std::vector<double> v_0(N, v_ini);

  for (int i = 1; i < N; i++) {
    double t = delta_t * i;
    v_0[i] = v[i - 1] + delta_t * dv_dt(v[i - 1], x[i - 1], (t - delta_t), m, k,
                                        c, f_0, w);
    x_0[i] = x[i - 1] + delta_t * v_0[i];
    v[i] = v[i - 1] +
           delta_t * dv_dt(v_0[i], x_0[i], (t - delta_t), m, k, c, f_0, w);
    x[i] = x[i - 1] + delta_t * v[i];
  }

  // Arquivo de saída
  if (save_to_file) {
    std::ofstream outfile(filename);

    if (outfile.is_open()) {
      outfile << "t,x,v" << std::endl;
      for (int i = 0; i < x.size(); i++) {
        outfile << i * delta_t << "," << x[i] << "," << v[i] << std::endl;
      }
      outfile.close();
      std::cout << "Simulação pelo métodod de Heun concluída com sucesso"
                << std::endl;
    } else {
      std::cout << "Não foi possível abrir o arquivo de saída" << std::endl;
    }
  } else {
    std::cout << "t,x,v" << std::endl;
    for (int i = 0; i < x.size(); i++) {
      std::cout << i * delta_t << "," << x[i] << "," << v[i] << std::endl;
    }
  }
}

int main() {
  double t_total = 20;   // s, tempo total de simulação
  double delta_t = 0.001; // s, passo de tempo
  double m = 30;          // kg, massa do sistema
  double c = 5;           // Ns/m, amortecimento viscoso
  double k = 1e3;         // N/m, rigidez do sistema
  double f_0 = 100;       // N, módulo do forçamento
  double w = 5;           // rad/s, frequência do forçamento

  integracao_euler(t_total, delta_t, m, c, k, f_0, w, 0.0, 0.0, true);
  integracao_heun(t_total, delta_t, m, c, k, f_0, w, 0.0, 0.0, true);

  return 0;
}
