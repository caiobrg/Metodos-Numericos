#include "solver_diferencas_finitas.hpp"
#include <cmath>
#include <cstdlib> // Necessário para o std::system
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

double f_lambda(double lb, double Bi) { return lb * std::tan(lb) - Bi; }

std::vector<double> calcular_autovalores(double Bi, int num_termos = 100) {
  std::vector<double> lambdas;
  double a;
  double b;
  double c;
  for (int i = 1; i <= num_termos; i++) {
    a = (i - 1) * M_PI + 1e-6;
    b = (i - 1) * M_PI + M_PI / 2.0 - 1e-6;
    // Bisseção simples
    for (int j = 0; j < 50; j++) {
      double c = (a + b) / 2.0;
      if (f_lambda(c, Bi) * f_lambda(a, Bi) < 0) {
        b = c;
      } else {
        a = c;
      }
    }
    lambdas.push_back((a + b) / 2.0);
  }
  return lambdas;
}

double theta_analitico(double x_star, double Fo,
                       const std::vector<double> &lbs) {
  // Retorna 1.0 imediatamente se o tempo (Fourier) for menor ou igual a zero
  if (Fo <= 0.0) {
    return 1.0;
  }

  double soma = 0.0;

  // Equivalente ao np.sum() e às operações vetorizadas do numpy
  for (double lambda : lbs) {
    double coef =
        (4.0 * std::sin(lambda)) / (2.0 * lambda + std::sin(2.0 * lambda));
    double termo =
        coef * std::cos(lambda * x_star) * std::exp(-(lambda * lambda) * Fo);
    soma += termo;
  }

  return soma;
}

void analise_convergencia(double L, double t_total, double t_inf,
                          double t_cond_inicial, double beta, double k,
                          double cp, double rho, double h1, double h2) {
  double geracao = 0.0; // Zero para o modelo analítico
  std::vector<int> ns_pontos = {4, 10, 40, 80, 160, 320, 500, 1000};
  std::vector<double> dts = {50.0, 10.0,  1.0,    0.5,
                             0.25, 0.125, 0.0625, 0.03125};

  // Parâmetros para a solução analítica
  double alpha = k / (cp * rho);
  double L_meio = L / 2.0;
  double Bi_ana = (h1 * L_meio) / k;
  double Fo_final = (alpha * t_total) / (L_meio * L_meio);

  // Pré-calcula os autovalores uma única vez
  std::vector<double> lambdas = calcular_autovalores(Bi_ana, 500);

  // ======================================================================
  // 1. ESTUDO DE CONVERGÊNCIA ESPACIAL (dt fixo, variando o número de nós)
  // ======================================================================
  std::cout << "\n--- ESTUDO DE CONVERGENCIA ESPACIAL (dt fixo = 0.001 s) ---"
            << std::endl;
  double dt_fixo = 0.001;

  // ABRE O ARQUIVO CSV ESPACIAL
  std::ofstream csv_esp("results/conv_espacial.csv");
  csv_esp << "N,R2\n";

  for (int n_pontos : ns_pontos) {
    std::string filename =
        "results/results_esp_N_" + std::to_string(n_pontos) + ".bin";
    double dx = L / (n_pontos - 1);

    // Executa o solver numérico
    solver_metodo_thomas(t_cond_inicial, t_inf, L, n_pontos, t_total, dt_fixo,
                         beta, geracao, 0.0, h1, h2, k, cp, rho, false, true,
                         filename.c_str());

    // Lê apenas o último frame (t_total) numérico gerado
    std::vector<double> T_num(n_pontos);
    std::ifstream arquivo_bin(filename, std::ios::binary);
    if (arquivo_bin) {
      arquivo_bin.seekg(-static_cast<long>(n_pontos * sizeof(double)),
                        std::ios::end);
      arquivo_bin.read(reinterpret_cast<char *>(T_num.data()),
                       n_pontos * sizeof(double));
      arquivo_bin.close();
    }

    // Calcula os valores analíticos exatos e as métricas do R^2
    std::vector<double> T_ana(n_pontos);
    double media_ana = 0.0;
    for (int i = 0; i < n_pontos; ++i) {
      double x_centro = (i * dx) - L_meio;
      double x_star = x_centro / L_meio;

      double theta = theta_analitico(x_star, Fo_final, lambdas);
      T_ana[i] = t_inf + theta * (t_cond_inicial - t_inf);
      media_ana += T_ana[i];
    }
    media_ana /= n_pontos;

    double ss_res = 0.0;
    double ss_tot = 0.0;
    for (int i = 0; i < n_pontos; ++i) {
      double res = T_num[i] - T_ana[i];
      ss_res += res * res;

      double tot = T_ana[i] - media_ana;
      ss_tot += tot * tot;
    }

    double r2 = (ss_tot > 0.0) ? (1.0 - (ss_res / ss_tot)) : 1.0;

    std::cout << "N = " << n_pontos << " \t| R^2 = " << r2 << std::endl;

    // SALVA NO CSV
    csv_esp << n_pontos << "," << std::fixed << std::setprecision(20) << r2
            << "\n";
  }
  csv_esp.close();

  // ======================================================================
  // 2. ESTUDO DE CONVERGÊNCIA TEMPORAL (Malha fixa, variando o dt)
  // ======================================================================
  std::cout << "\n--- ESTUDO DE CONVERGENCIA TEMPORAL (N fixo = 1000 nós) ---"
            << std::endl;
  int n_fixo = 1000;
  double dx_fixo = L / (n_fixo - 1);

  // ABRE O ARQUIVO CSV TEMPORAL
  std::ofstream csv_temp("results/conv_temporal.csv");
  csv_temp << "dt,R2\n";

  // A malha espacial é fixa, logo a analítica final é idêntica para todos os
  // dt's. Podemos pré-calcular para economizar processamento.
  std::vector<double> T_ana_temp(n_fixo);
  double media_ana_temp = 0.0;
  for (int i = 0; i < n_fixo; ++i) {
    double x_centro = (i * dx_fixo) - L_meio;
    double x_star = x_centro / L_meio;

    double theta = theta_analitico(x_star, Fo_final, lambdas);
    T_ana_temp[i] = t_inf + theta * (t_cond_inicial - t_inf);
    media_ana_temp += T_ana_temp[i];
  }
  media_ana_temp /= n_fixo;

  for (double dt_teste : dts) {
    std::ostringstream oss;
    oss << "results/results_temp_dt_" << dt_teste << ".bin";
    std::string filename = oss.str();

    // Executa o solver numérico
    solver_metodo_thomas(t_cond_inicial, t_inf, L, n_fixo, t_total, dt_teste,
                         beta, geracao, 0.0, h1, h2, k, cp, rho, false, true,
                         filename.c_str());

    // Lê o último frame numérico
    std::vector<double> T_num(n_fixo);
    std::ifstream arquivo_bin(filename, std::ios::binary);
    if (arquivo_bin) {
      arquivo_bin.seekg(-static_cast<long>(n_fixo * sizeof(double)),
                        std::ios::end);
      arquivo_bin.read(reinterpret_cast<char *>(T_num.data()),
                       n_fixo * sizeof(double));
      arquivo_bin.close();
    }

    // Calcula R^2
    double ss_res = 0.0;
    double ss_tot = 0.0;
    for (int i = 0; i < n_fixo; ++i) {
      double res = T_num[i] - T_ana_temp[i];
      ss_res += res * res;

      double tot = T_ana_temp[i] - media_ana_temp;
      ss_tot += tot * tot;
    }

    double r2 = (ss_tot > 0.0) ? (1.0 - (ss_res / ss_tot)) : 1.0;
    csv_temp << dt_teste << "," << std::fixed << std::setprecision(40) << r2
             << "\n";
    std::cout << "dt = " << dt_teste << " s \t| R^2 = " << r2 << std::endl;
  }
  csv_temp.close();

  // ======================================================================
  // 3. CHAMADA DO PYTHON PARA PLOTAR
  // ======================================================================
  std::cout << "\nGerando graficos de convergencia via Python..." << std::endl;
  std::ostringstream comando_python;
  comando_python << "python3 plot_convergencia.py "
                 << "--arq_esp results/conv_espacial.csv "
                 << "--arq_temp results/conv_temporal.csv "
                 << "--L " << L;

  std::system(comando_python.str().c_str());
}

void validacao_animada(double L, int n_pontos, double t_total, double dt,
                       double t_inf, double t_cond_inicial, double beta,
                       double k, double cp, double rho, double h1, double h2) {
  std::cout << "--- Iniciando Simulação de Validação ---" << std::endl;
  double geracao = 0.0; // Zero para o modelo analítico

  std::string filename = "results/results_validacao_animada.bin";
  // --- Executar o Solver Numérico ---
  solver_metodo_thomas(t_cond_inicial, t_inf, L, n_pontos, t_total, dt, beta,
                       geracao, 0.0, h1, h2, k, cp, rho, false, true,
                       filename.c_str());

  std::ostringstream comando_python;
  comando_python << "python3 plot_validacao.py " << "--arquivo " << filename
                 << " --tempo_total " << t_total << " --dt " << dt
                 << " --pontos " << n_pontos << " --L " << L << " --k " << k
                 << " --cp " << cp << " --rho " << rho << " --h1 " << h1
                 << " --t_inf " << t_inf << " --t_cond_inicial "
                 << t_cond_inicial;

  std::cout << "Executando animação de validação..." << std::endl;
  std::system(comando_python.str().c_str());
}

void simular_e_animar_heatmap(double t_cond_ini, double t_inf, double L,
                              int n_pontos, double t_total, double dt,
                              double beta, double geracao, double k, double cp,
                              double rho, double h1, double h2,
                              std::string titulo) {
  std::cout << "\n--- Iniciando Simulação para Heatmap: " << titulo << " ---"
            << std::endl;

  std::string filename = "results/results_heatmap.bin";
  solver_metodo_thomas(t_cond_ini, t_inf, L, n_pontos, t_total, dt, beta,
                       geracao, 0.0, h1, h2, k, cp, rho, false, true,
                       filename.c_str());

  std::ostringstream comando;
  comando << "python3 heatmap_animation.py "
          << "--arquivo " << filename << " --tempo_total " << t_total
          << " --dt " << dt << " --pontos " << n_pontos << " --titulo \""
          << titulo << "\"";

  std::cout << "Executando: " << comando.str() << std::endl;
  std::system(comando.str().c_str());
}

int main() {
  // --- Parâmetros Físicos e Numéricos ---
  double L = 0.01;
  int n_pontos = 100;
  double t_total = 100.0;
  double dt = 0.125;
  double t_inf = 300;
  double t_cond_inicial = 300;
  double beta = 0.5;

  double k = 4.0;
  double cp = 490.0;
  double rho = 10500.0;
  double h1 = 10000.0;
  double h2 = 30000.0;

  // Chamada da análise de convergência de malha
  analise_convergencia(L, 100, 20, 300, beta, k, cp, rho, 40, 40);

  // Chamada da Validação (Gráfico de Curvas Animado)
  validacao_animada(L, n_pontos, t_total, dt, t_inf, 1000, beta, k, cp, rho, h1,
                    h1);

  // Chamada do Heatmap Animado com Geração Interna
  double geracao = 3.0e8; // W/m^3
  simular_e_animar_heatmap(t_cond_inicial, t_inf, L, n_pontos, t_total, dt,
                           beta, geracao, k, cp, rho, h1, h2,
                           "Heatmap - Resfriamento com Geracao Interna");

  return 0;
}
