#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <ratio>
#include <string>
#include <vector>

class configSimulacao {
private:
  double L;    // Comprimento da aleta [m]
  double H;    // Altura da aleta [m]
  int Nx;      // Número de nós na direção x
  int Ny;      // Número de nós na direção y
  double Tb;   // Temperatura na base [K]
  double Tinf; // Temperatura nas demais parede [K]
  double tol;  // Tolerância de convergência
  int maxIt;   // Número máximo de iterações

  double dx;
  double dy;
  double beta;

  // Método auxiliar privado para captura de input (sobrecarga para double)
  double capturarInput(double valorPadrao, const std::string &nome) {
    std::string input;

    std::cout << "Digite o valor de " << nome << " [default: " << valorPadrao
              << "]: ";
    std::getline(std::cin, input);

    if (!input.empty()) {
      return std::stod(input);
    }
    return valorPadrao;
  }

  // Método auxiliar privado para captura de input (sobrecarga para int)
  int capturarInput(int valorPadrao, const std::string &nome) {
    std::string input;
    std::cout << "Digite o valor de " << nome << " [default: " << valorPadrao
              << "]: ";
    std::getline(std::cin, input);

    if (!input.empty()) {
      return std::stoi(input); // Usa stoi especificamente para inteiros
    }
    return valorPadrao;
  }

public:
  // Construção da classe
  configSimulacao(int Nx = 41, int Ny = 41, double L = 0.5, double H = 1,
                  double Tb = 300, double Tinf = 273, double tol = 1e-4,
                  int maxIt = 200000)
      : L(L), H(H), Nx(Nx), Ny(Ny), Tb(Tb), Tinf(Tinf), tol(tol), maxIt(maxIt) {

    dx = this->L / (Nx - 1);
    dy = this->H / (Ny - 1);
    beta = std::pow(dx / dy, 2);
  }
  // Prompta o usuário no terminal sobre os valores desejados
  void configurarTerminal() {
    L = capturarInput(L, "Comprimento da aleta (L) [m]");
    H = capturarInput(H, "Altrua da aleta (H) [m]");
    Nx = capturarInput(Nx, "Número de pontos em L (Nx)");
    Ny = capturarInput(Ny, "Número de pontos em H (Ny)");
    Tb = capturarInput(Tb, "Temperatura da base da aleta (Tb) [K]");
    Tinf = capturarInput(Tinf,
                         "Temperatura das demais paredes da aleta (T_inf) [K]");
    tol = capturarInput(tol, "Tolerância para o método iterativo (tol)");
    maxIt = capturarInput(maxIt, "Número máximo de iterações (maxIt)");
  }

  // Metodos getters, retornam os valores em modo somente leitura
  double get_L() const { return L; }
  double get_H() const { return H; }
  int get_Nx() const { return Nx; }
  int get_Ny() const { return Ny; }
  double get_T_b() const { return Tb; }
  double get_T_inf() const { return Tinf; }
  double get_tol() const { return tol; }
  double get_max_It() const { return maxIt; }
  double get_dx() const { return dx; }
  double get_dy() const { return dy; }
  double get_beta() const { return beta; }
};

class Node {
private:
  // Variáveis imutáveis do nó
  size_t id; // Id único para identificação do nó
  double x;  // Posição do nó em x
  double y;  // Posição do nó em y
public:
  double T; // Temperatura do nó, pode ser alterada durante a execução do nó
  Node(int id, double x, double y, double T) : id(id), x(x), y(y), T(T) {}

  // Metodos getters, retornam os valores em modo somente leitura
  double get_id() const { return id; }
  double get_x() const { return x; }
  double get_y() const { return y; }
};

class Mesh {
  configSimulacao config;  // Configuração a ser passada na criação do objeto
  std::vector<Node> nodes; // Vetor para guardar todos os nós

public:
  Mesh(configSimulacao configSimulacao) : config(configSimulacao) {}
  void buildMesh() {
    double Tinf = config.get_T_inf();
    double Tb = config.get_T_b();

    double Nx = config.get_Nx();
    double Ny = config.get_Ny();

    double dx = config.get_dx();
    double dy = config.get_dy();

    // Constróid a malha de fato
    size_t Id = 0;
    for (double pos_y = 0; pos_y < Ny; pos_y++) {
      double y = pos_y * dy;
      for (double pos_x = 0; pos_x < Nx; pos_x++) {
        double x = pos_x * dx;
        if (pos_y == 0) {
          nodes.push_back(Node(Id, x, y, Tb));
        } else {
          nodes.push_back(Node(Id, x, y, Tinf));
        }
        Id++;
      }
    }
  }

  std::vector<Node> get_Nodes() const { return nodes; }

  void save_mesh(std::string filename, bool verbose = false) {
    std::ofstream meshfile(filename);

    // Printa o cabeçalho no arquivo de saída
    meshfile << "id" << "," << "x" << ","
             << "y" << "," << "T" << std::endl;

    // Percorre a malha nó por nó e salva um por linha no arquivo de sáida
    for (const Node &node : nodes) {
      meshfile << std::setprecision(12) << node.get_id() << "," << node.get_x()
               << "," << node.get_y() << "," << node.T << std::endl;
    }
    // fecha o arquivo de saída e avisa o usuário do fim do salvamento
    meshfile.close();
    std::cout << "Malha salva no arquivo " << filename << std::endl;

    if (verbose) {
      std::cout << "id" << std::setw(15) << "x" << std::setw(15) << "y"
                << std::setw(15) << "T" << std::endl;

      for (const Node &node : nodes) {
        std::cout << std::setprecision(6) << node.get_id() << std::setw(15)
                  << node.get_x() << std::setw(15) << node.get_y()
                  << std::setw(15) << node.T << std::endl;
      }
    }
  }

  // Mapeia a posição 2D (i.j) para um array 1D (j * Nx + i)
  Node &get_Node(int i, int j) {
    int Nx = config.get_Nx();
    return nodes[j * Nx + i];
  }
};

class sistema {
private:
  Mesh &mesh;
  configSimulacao &config;

  double beta = config.get_beta();

  double dar_passo() {
    double resMax = 0;
    for (size_t i = 1; i < config.get_Nx() - 1; i++) {
      for (size_t j = 1; j < config.get_Ny() - 1; j++) {
        Node &node = mesh.get_Node(i, j);
        Node &nodeS = mesh.get_Node(i, j - 1);
        Node &nodeN = mesh.get_Node(i, j + 1);
        Node &nodeE = mesh.get_Node(i + 1, j);
        Node &nodeW = mesh.get_Node(i - 1, j);
        double T_old = node.T;
        node.T =
            (nodeW.T + nodeE.T + beta * (nodeN.T + nodeS.T)) / (2 * (1 + beta));

        double res = std::abs(node.T - T_old);

        if (res > resMax)
          resMax = res;
      }
    }
    return resMax;
  }

public:
  sistema(Mesh &mesh, configSimulacao &config) : mesh(mesh), config(config) {}

  void solve(bool verbose = true) {
    double res = 1;
    size_t it = 0;
    auto t1 = std::chrono::high_resolution_clock::now();
    while (res >= config.get_tol() && it <= config.get_max_It()) {
      res = dar_passo();
      it++;
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::milli> time_ms = t2 - t1;

    if (verbose) {
      std::cout << "Número de iterações para a convergência: " << it
                << std::endl;
      std::cout << "Resíduo máximo da última iteração: " << std::scientific
                << res << std::endl;
      std::cout << "Tempo de execução do loop do solver: " << time_ms.count()
                << std::setprecision(4) << " ms" << std::endl;
    }
  }
};

int main() {
  configSimulacao arrayConfig[]{configSimulacao{11, 11}, configSimulacao{21, 21}, configSimulacao{41, 41},
                                 configSimulacao{81, 81}, configSimulacao{500, 500}};
  for (configSimulacao &configuracao : arrayConfig) {
    Mesh malha(configuracao);

    malha.buildMesh();
    malha.save_mesh(std::string("results/malha_") +
                    std::to_string((configuracao.get_Nx())) + "x" +
                    std::to_string(configuracao.get_Ny()) + ".csv");

    sistema sistema(malha, configuracao);
    sistema.solve();
    malha.save_mesh(std::string("results/malha_resolvida_") +
                    std::to_string(configuracao.get_Nx()) + "x" +
                    std::to_string(configuracao.get_Ny()) + ".csv");
    std::cout << std::string(70, '-') << std::endl;
  }
}
