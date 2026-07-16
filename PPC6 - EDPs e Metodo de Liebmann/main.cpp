#include "classes.hpp"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>

int main() {
  // Instancia a configuração e a malha usando os valores defaults
  ConfigSimulacao config(21, 21);
  Mesh malha{config};
  malha.build_mesh();

  // Configura as condições de contorno
  System sistema(malha);

  // Demais arestas: Robin (Tinf, h)
  sistema.set_bc(malha.north()).robin(config.get_T_inf(), config.get_h());
  sistema.set_bc(malha.south()).robin(config.get_T_inf(), config.get_h());
  sistema.set_bc(malha.east()).robin(config.get_T_inf(), config.get_h());

  // Aresta Oeste (x=0) é a base: Dirichlet (Tb)
  sistema.set_bc(malha.west()).dirichlet(config.get_T_b());

  // Instancia o solver
  Solver solver(sistema);

  // 1. Eliminação de Gauss
  std::cout << "--- 1. Eliminacao de Gauss ---\n";
  Solver::SolverResult res_gauss = solver.eliminacao_gaussiana();
  std::cout << "Tempo de execucao: " << res_gauss.time_ms << " ms\n";
  sistema.salvar_resultados("resultados/resultados_gauss.csv");

  sistema.reset_temperatures();

  // 2. Liebmann sem relaxamento (omega = 1.0)
  std::cout << "\n--- 2. Metodo de Liebmann (sem relaxamento) ---\n";
  Solver::SolverResult res_lieb = solver.liebmann(1.0);
  std::cout << "Iteracoes: " << res_lieb.iterations << "\n";
  std::cout << "Erro final: " << res_lieb.erro_max << "\n";
  std::cout << "Tempo de execucao: " << res_lieb.time_ms << " ms\n";
  sistema.salvar_resultados("resultados/resultados_liebmann.csv");

  sistema.reset_temperatures();

  // 3. Liebmann com relaxamento (omega = 1.5)
  std::cout << "\n--- 3. Metodo de Liebmann (com relaxamento, w = 1.5) ---\n";
  Solver::SolverResult res_lieb_relax = solver.liebmann(1.5);
  std::cout << "Iteracoes: " << res_lieb_relax.iterations << "\n";
  std::cout << "Erro final: " << res_lieb_relax.erro_max << "\n";
  std::cout << "Tempo de execucao: " << res_lieb_relax.time_ms << " ms\n";
  sistema.salvar_resultados("resultados/resultados_liebmann_relax.csv");

  // -------------------------------------------------------------
  // ESTUDO 1: Efeito do relaxamento na malha 81x81
  // -------------------------------------------------------------
  std::cout << "\n===================================================\n";
  std::cout << "INICIANDO ESTUDO DO FATOR DE RELAXACAO OMEGA (81x81)\n";
  std::cout << "===================================================\n";

  ConfigSimulacao config_relax(81, 81);
  Mesh malha_relax(config_relax);
  malha_relax.build_mesh();

  std::ofstream file_relax("resultados/estudo_relaxacao.csv");
  file_relax << "omega,iteracoes,tempo_ms\n";

  for (double w = 1.0; w < 1.95; w += 0.05) {
    System sistema_relax(malha_relax);

    sistema_relax.set_bc(malha_relax.north())
        .robin(config_relax.get_T_inf(), config_relax.get_h());
    sistema_relax.set_bc(malha_relax.south())
        .robin(config_relax.get_T_inf(), config_relax.get_h());
    sistema_relax.set_bc(malha_relax.east())
        .robin(config_relax.get_T_inf(), config_relax.get_h());

    sistema_relax.set_bc(malha_relax.west()).dirichlet(config_relax.get_T_b());

    Solver solver_relax(sistema_relax);
    Solver::SolverResult res = solver_relax.liebmann(w);

    std::cout << "omega = " << std::fixed << std::setprecision(2) << w
              << " | Iteracoes: " << res.iterations
              << " | Tempo: " << std::fixed << std::setprecision(4)
              << res.time_ms << " ms\n";

    file_relax << std::fixed << std::setprecision(2) << w << ","
               << res.iterations << "," << std::fixed << std::setprecision(4)
               << res.time_ms << "\n";
  }
  file_relax.close();
  std::cout << "Resultados salvos em 'resultados/estudo_relaxacao.csv'\n";

  // -------------------------------------------------------------
  // ESTUDO 2: Convergência de malha
  // -------------------------------------------------------------
  std::cout << "\n===================================================\n";
  std::cout << "INICIANDO ESTUDO DE REFINAMENTO DE MALHA\n";
  std::cout << "===================================================\n";

  std::vector<int> tamanhos = {11, 21, 41, 81};

  for (int N : tamanhos) {
    ConfigSimulacao config_conv(N, N);
    Mesh malha_conv(config_conv);
    malha_conv.build_mesh();

    System sistema_conv(malha_conv);

    sistema_conv.set_bc(malha_conv.north())
        .robin(config_conv.get_T_inf(), config_conv.get_h());
    sistema_conv.set_bc(malha_conv.south())
        .robin(config_conv.get_T_inf(), config_conv.get_h());
    sistema_conv.set_bc(malha_conv.east())
        .robin(config_conv.get_T_inf(), config_conv.get_h());
    sistema_conv.set_bc(malha_conv.west()).dirichlet(config_conv.get_T_b());

    Solver solver_conv(sistema_conv);

    // Liebmann com relaxação 1.5
    Solver::SolverResult res_conv = solver_conv.liebmann(1.5);

    std::cout << "Malha " << N << "x" << N
              << " | Iteracoes: " << res_conv.iterations
              << " | Tempo: " << std::fixed << std::setprecision(4)
              << res_conv.time_ms << " ms\n";

    std::string filename = "resultados/temperatura_" + std::to_string(N) + "x" +
                           std::to_string(N) + ".csv";
    sistema_conv.salvar_resultados(filename);
  }

  return 0;
}