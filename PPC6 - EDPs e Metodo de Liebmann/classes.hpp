#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

// Classe que armazena as constantes que serão usadas para toda a simulação
class ConfigSimulacao {
private:
  double L; // Comprimento da aleta [m]
  double H; // Altura da aleta [m]
  int Nx;   // Número de nós na direção x
  int Ny;   // Número de nós na direção y

  double TInf; // Temperatura do escoamento [K] padrão,
               // não é usada por default no código

  double Tb; // Temperatura da base aquecida [K] padrão,
             // não é usada por default no código

  double k; // Condutividade térmica [W/(m * K)]
  double h; // Coeficiente convectivo padrão [W/(m^2 * K)]

  double tol; // Tolerância de convergência
  int maxIt;  // Número máximo de iterações

  // Constantes geométricas
  double dx;   // Espaçamento entre os nós na direção x [m]
  double dy;   // Espaçamento entre os nós na direção y [m]
  double beta; // Razão quadrática dos espaçamentos

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
  // Construtor da classe com valores padrão
  ConfigSimulacao(int Nx = 11, int Ny = 11, double L = 0.5, double H = 1,
                  double k = 200, double h = 100, double Tb = 425,
                  double Tinf = 273, double tol = 1e-4, int maxIt = 200000)
      : L(L), H(H), Nx(Nx), Ny(Ny), Tb(Tb), TInf(Tinf), k(k), h(h), tol(tol),
        maxIt(maxIt) {

    // Calcula as constantes geométricas, this só especifica que o H e L são oo
    // passados para a função
    dx = this->L / (Nx - 1);
    dy = this->H / (Ny - 1);
    beta = std::pow(dx / dy, 2);
  }
  // Prompta o usuário no terminal sobre os valores desejados, caso ele não
  // insira nenhum é usado o valor default
  void configurarTerminal() {
    L = capturarInput(L, "Comprimento da aleta (L) [m]");
    H = capturarInput(H, "Altura da aleta (H) [m]");
    Nx = capturarInput(Nx, "Número de pontos em L (Nx)");
    Ny = capturarInput(Ny, "Número de pontos em H (Ny)");
    tol = capturarInput(tol, "Tolerância para o método iterativo (tol)");
    maxIt = capturarInput(maxIt, "Número máximo de iterações (maxIt)");
    k = capturarInput(k, "Condutividade térmica do material [W/(m * K)]");
    h = capturarInput(h, "Coeficiente convectivo padrão [W/(m^2 * K)]");
    Tb = capturarInput(Tb, "Temperatura padrão da base aquecida [K]");
    TInf = capturarInput(TInf, "Temperatura padrão do escoamento [K]");

    // Recalcula as constantes geométricas caso ocorra configuração pelo
    // terminal
    dx = this->L / (Nx - 1);
    dy = this->H / (Ny - 1);
    beta = std::pow(dx / dy, 2);
  }

  // métodos getters, retornam os valores em modo somente leitura
  double get_L() const { return L; }
  double get_H() const { return H; }
  int get_Nx() const { return Nx; }
  int get_Ny() const { return Ny; }
  double get_k() const { return k; }
  double get_h() const { return h; }
  double get_T_b() const { return Tb; }
  double get_T_inf() const { return TInf; }
  double get_tol() const { return tol; }
  int get_max_It() const { return maxIt; }
  double get_dx() const { return dx; }
  double get_dy() const { return dy; }
  double get_beta() const { return beta; }
};

// Classe representativa de um nó, guarda somente as características geométricas
class Node {
private:
  // Variáveis imutáveis do nó
  size_t id; // Id único para identificação do nó
  double x;  // Posição do nó em x
  double y;  // Posição do nó em y
public:
  Node(int id, double x, double y) : id(id), x(x), y(y) {}

  // métodos getters, retornam os valores em modo somente leitura
  double get_id() const { return id; }
  double get_x() const { return x; }
  double get_y() const { return y; }
};

// Declaração adiantada apenas para que possa se usar friend class dentro de
// Mesh
class System;

// Classe para a construção da malha, ainda não possui as condições de
// contorno e nem os valores de temperatura
class Mesh {
  const ConfigSimulacao &config; // Referência constante à configuração original
                                 // (impede que se altere a config original)

  std::vector<Node> meshnodes; // Vetor para guardar todos os nós

  // Permite que um objeto da classe System acesse as variáveis privadas
  friend class System;

public:
  // Construtor da classe
  Mesh(const ConfigSimulacao &configSimulacao) : config(configSimulacao) {}

  // Overload do operador (), permite extrair um nó da malha usando Mesh(i,j)
  Node &operator()(int i, int j) { return meshnodes[i * config.get_Nx() + j]; }

  // Overload do operador (), permite extrair um nó da malha usando Mesh(i,j),
  // versão para funções que exijam read-only
  const Node &operator()(int i, int j) const {
    return meshnodes[i * config.get_Nx() + j];
  }

  // Função para calcular a malha
  void build_mesh() {

    // Extração das variáveis da configuração para melhorar a sintaxe
    double Nx = config.get_Nx();
    double Ny = config.get_Ny();

    double dx = config.get_dx();
    double dy = config.get_dy();

    // Constrói a malha com as temperaturas inicializadas em 0 K
    size_t Id = 0;
    // Varredura em y
    for (double pos_y = 0; pos_y < Ny; pos_y++) {
      double y = pos_y * dy; // Posição y em m para coordenada do nó
      // Varredura em x para
      for (double pos_x = 0; pos_x < Nx; pos_x++) {
        double x = pos_x * dx; // Posição x em m para coordenada do nó
        meshnodes.push_back(
            Node(Id, x, y)); // Adiciona o nó ao vetor de coordenada
        Id++;
      }
    }
  }
  const int get_total_size() const { return meshnodes.size(); }

  // Mesh.south() retorna todos os ids dos nós da fronteira Sul da malha
  std::vector<int> south() {
    int Nx = config.get_Nx();
    // South nodes é um vetor de ids (inteiros)
    std::vector<int> southNodes;

    // Varre cada coluna da primeira linha e adiciona ao vetor de referências
    for (int j = 0; j < Nx; j++) {
      // (*this) referencia o objeto que chamou o método, logo Mesh(0,j)
      southNodes.push_back((*this)(0, j).get_id());
    }
    return southNodes;
  }
  // Mesh.north() retorna todos os ids dos nós da fronteira Norte da malha, para
  // comentários mais completos se referir ao método south()
  std::vector<int> north() {
    int Ny = config.get_Ny();
    int Nx = config.get_Nx();
    std::vector<int> northNodes;
    for (int j = 0; j < Nx; j++) {
      northNodes.push_back((*this)(Ny - 1, j).get_id());
    }
    return northNodes;
  }
  // Mesh.west() retorna todos os ids dos nós da fronteira Oeste da malha, para
  // comentários mais completos se referir ao método south()
  std::vector<int> west() {
    int Ny = config.get_Ny();
    int Nx = config.get_Nx();
    std::vector<int> westNodes;
    for (int i = 0; i < Ny; i++) {
      westNodes.push_back((*this)(i, 0).get_id());
    }
    return westNodes;
  }
  // Mesh.east() retorna todos os nós da fronteira Leste da malha, para
  // comentários mais completos se referir ao método south()
  std::vector<int> east() {
    int Ny = config.get_Ny();
    int Nx = config.get_Nx();
    std::vector<int> eastNodes;
    for (int i = 0; i < Ny; i++) {
      eastNodes.push_back((*this)(i, Nx - 1).get_id());
    }
    return eastNodes;
  }

  // Função para salvar os nós da malha, salva um nó por linha com suas
  // informações separadas por ","
  void save_mesh(std::string filename, bool verbose = false) {
    std::ofstream meshfile(filename);

    // Printa o cabeçalho no arquivo de saída
    meshfile << "id" << "," << "x" << ","
             << "y" << "," << std::endl;

    // Percorre a malha nó por nó e salva um por linha no arquivo de saída
    for (const Node &node : meshnodes) {
      meshfile << std::setprecision(12) << node.get_id() << "," << node.get_x()
               << "," << node.get_y() << std::endl;
    }
    // fecha o arquivo de saída e avisa o usuário do fim do salvamento
    meshfile.close();
    std::cout << "Malha salva no arquivo " << filename << std::endl;

    if (verbose) {
      std::cout << "id" << std::setw(15) << "x" << std::setw(15) << "y"
                << std::endl;

      for (const Node &node : meshnodes) {
        std::cout << std::setprecision(6) << node.get_id() << std::setw(15)
                  << node.get_x() << std::setw(15) << node.get_y() << std::endl;
      }
    }
  }
};

// Declaração adiantada apenas para que possa se usar friend class dentro de
// System
class Solver;

// Classe do sistema físico, guarda temperaturas, condições de contorno. Ainda
// não soluciona
class System {
private:
  const Mesh &malha; // Referência constante à malha (impede que se altere a
                     // malha original)

  const ConfigSimulacao &config; // Referência constante à configuração
                                 // (será extraída da malha)

  // Classe que guarda define os tipos de condição, para evitar o uso de strings
  enum class BcTypes { None, Dirichlet, Robin };

  // Vetores de estado físico indexados pelo ID do nó
  std::vector<double> temperatures;
  std::vector<BcTypes>
      bcTypes; // BcTypes::Dirichlet, BcTypes::Robin, BcTypes::None
  std::vector<double> bcValues; // Guarda Tb ou T_inf
  std::vector<double> hCoeffs;  // Guarda o h para Robin

  // Classe aninhada que recebe os IDs alvo e aplica a condição de contorno.
  // Permite que seja usado System.set_bc(IDs).robin(T_inf, h)
  class BoundaryProxy {
  private:
    std::vector<int> target_ids; // Cópia do vetor de ids a ser alterado
    System &sistema;             // Referência

  public:
    // Construtor para o caso em que é passado um vetor de ids (ex:
    // mesh.south)
    BoundaryProxy(const std::vector<int> &ids, System &sistema)
        : target_ids(ids), sistema(sistema) {}

    BoundaryProxy(int id, System &sistema) : target_ids{id}, sistema(sistema) {}

    // helper function para avisar caso haja conflito em algum nó
    void _avisar_conflito(int id) {
      // Código ANSI \033[1;33m deixa o texto amarelo no
      // terminal, e \033[0m reseta a cor
      std::cerr << "\033[1;33m[AVISO] Conflito no nó ID " << id
                << ": mantendo a última condição especificada\033[0m\n";
    }

    // Sobrecarga 1: Dirichlet com valor único uniforme
    void dirichlet(double Tb) {

      // Itera sobre cada um dos IDs alvo
      for (int id : target_ids) {

        // Verifica se o nó já havia sido definido como robin
        if (sistema.bcTypes[id] == BcTypes::Robin) {
          _avisar_conflito(id);
        }

        // Define o tipo de condição de contorno e a temperatura fixa
        sistema.bcTypes[id] = BcTypes::Dirichlet;
        sistema.bcValues[id] = Tb;
      }
    }

    // Sobrecarga 2: Dirichlet com um vetor de valores
    void dirichlet(const std::vector<double> &Tb_vector) {

      // Verifica se o vetor tem o mesmo tamanho do número de IDs alvo, evita
      // erro de out of bounds
      if (Tb_vector.size() != target_ids.size()) {
        throw std::runtime_error(
            "O tamanho do vetor de temperaturas nao condiz "
            "com o numero de nos na fronteira!");
      }
      for (size_t i = 0; i < target_ids.size(); ++i) {
        int id = target_ids[i];

        // Verifica se o nó já havia sido definido como robin
        if (sistema.bcTypes[id] == BcTypes::Robin) {
          _avisar_conflito(id);
        }

        // Define o tipo de condição de contorno e a temperatura fixa
        sistema.bcTypes[id] = BcTypes::Dirichlet;
        sistema.bcValues[id] = Tb_vector[i];
      }
    }

    // Sobrecarga 3: Robin com valor único uniforme uniforme
    void robin(double T_inf, double h) {
      for (int id : target_ids) {

        // Verifica se o nó já tinha sido definido como Dirichlet anteriormente
        if (sistema.bcTypes[id] == BcTypes::Dirichlet) {
          _avisar_conflito(id);
        }

        // Define o tipo de condição de contorno,a temperatura do escoamento e o
        // coeficiente convectivo
        sistema.bcTypes[id] = BcTypes::Robin;
        sistema.bcValues[id] = T_inf;
        sistema.hCoeffs[id] = h;
      }
    }

    // Sobrecarga 4: Robin com vetores de valores
    void robin(const std::vector<double> &Tinf_vector,
               const std::vector<double> &h_vector) {

      // Verifica se os vetores tem o mesmo tamanho do número de IDs alvo, evita
      // erro de out of bounds
      if (Tinf_vector.size() != target_ids.size()) {
        throw std::runtime_error(
            "O tamanho do vetor de temperaturas nao condiz "
            "com o numero de nos na fronteira!");
      } else if (h_vector.size() != target_ids.size()) {
        throw std::runtime_error(
            "O tamanho do vetor de coeficientes convectivos nao condiz "
            "com o numero de nos na fronteira!");
      }

      for (size_t i = 0; i < target_ids.size(); ++i) {
        int id = target_ids[i];
        // Verifica se o nó já tinha sido definido como Dirichlet anteriormente
        if (sistema.bcTypes[id] == BcTypes::Dirichlet) {
          _avisar_conflito(id);
        }

        // Define o tipo de condição de contorno,a temperatura do escoamento e o
        // coeficiente convectivo
        sistema.bcTypes[id] = BcTypes::Robin;
        sistema.bcValues[id] = Tinf_vector[i];
        sistema.hCoeffs[id] = h_vector[i];
      }
    }
  };

  // Permite que um objeto da classe solver acesse as variáveis privadas
  friend class Solver;

public:
  // Construtor da classe, inicializa os vetores com seu tamanho alvo e valores
  // nulos
  System(const Mesh &malha)
      : malha(malha), config(malha.config),
        temperatures(malha.get_total_size(), 0.0),
        bcTypes(malha.get_total_size(), BcTypes::None),
        bcValues(malha.get_total_size(), 0.0),
        hCoeffs(malha.get_total_size(), 0.0) {}

  // Função para adicionar uma condição de contorno ao sistema. Overload para
  // int (um único ID)
  BoundaryProxy set_bc(int node_id) { return BoundaryProxy(node_id, (*this)); }

  // Overload para um vetor de int (múltiplos IDs)
  BoundaryProxy set_bc(std::vector<int> node_id) {
    return BoundaryProxy(node_id, (*this));
  }

  // Reseta o campo de temperatura para o valor inicial
  void reset_temperatures(double initial_T = 0.0) {
    std::fill(temperatures.begin(), temperatures.end(), initial_T);
  }

  // Salva os resultados (id, x, y, T) em um arquivo CSV
  void salvar_resultados(std::string filename) {
    std::ofstream meshfile(filename);
    meshfile << "id,x,y,T\n";
    // Como System é friend de Mesh, podemos acessar meshnodes diretamente
    for (size_t i = 0; i < malha.get_total_size(); i++) {
      const Node &node = malha(i / config.get_Nx(), i % config.get_Nx());
      meshfile << std::setprecision(12) << node.get_id() << "," << node.get_x()
               << "," << node.get_y() << "," << temperatures[node.get_id()]
               << "\n";
    }
    meshfile.close();
    std::cout << "Resultados salvos no arquivo " << filename << std::endl;
  }
};

class Solver {
private:
  System &sistema;

  // Struct para guardar uma linha da matriz
  struct Equacao {
    double diag;
    std::vector<std::pair<int, double>>
        vizinhos; // (id do vizinho, coeficiente)
    double rhs;
  };

  Equacao monta_equacao(int i) {
    // Características da malha
    int Nx = sistema.config.get_Nx();
    int Ny = sistema.config.get_Ny();
    int ix = i % Nx, iy = i / Nx;

    // Condição de contorno do nó analisado
    System::BcTypes bc = sistema.bcTypes[i];

    // Objeto à ser retornado
    Equacao eq;

    // Se  a condição de contorno é dirichlet
    if (bc == System::BcTypes::Dirichlet) {
      // T_ij = Tb
      eq.diag = 1.0;
      eq.rhs = sistema.bcValues[i];
      return eq; // sem vizinhos
    }

    // Se a condição de contorno é robin
    if (bc == System::BcTypes::Robin) {
      // Características físicas do sistema
      double k = sistema.config.get_k(); // Condutividade do material
      double h = sistema.hCoeffs[i];     // Coef convectivo naquele nó
      double dx = sistema.config.get_dx(),
             dy = sistema.config.get_dy(); // Espaçamentos

      // Agrupamentos para facilitar a notação
      double theta = (2 * dx * h) / k, omega = (2 * dy * h) / k;
      double beta = sistema.config.get_beta();
      double sb = std::pow(beta, 0.5), isb = std::pow(beta, -0.5);

      // Verifica se o nó se encontra em algum dos cantos
      bool is_corner = (ix == 0 || ix == Nx - 1) && (iy == 0 || iy == Ny - 1);

      if (is_corner) {
        // Outro agrupamento para facilitar a notação
        double lambda = (theta + omega) / 2.0;

        eq.diag = -(sb + isb + lambda);

        // Verifica qual a posição do nó (importante para definir quais os
        // vizinhos)
        if (ix != 0 && iy != 0) { // nó NE
          eq.vizinhos = {{i - 1, isb}, {i - Nx, sb}};
        } else if (ix == 0 && iy == 0) { // nó SW
          eq.vizinhos = {{i + 1, isb}, {i + Nx, sb}};
        } else if (ix != 0 && iy == 0) { // nó SE
          eq.vizinhos = {{i - 1, isb}, {i + Nx, sb}};
        } else { // nó NW
          eq.vizinhos = {{i + 1, isb}, {i - Nx, sb}};
        }

        eq.rhs = -lambda * sistema.bcValues[i];

      } else if (iy == 0 || iy == Ny - 1) { // Aresta inferior/superior
        // Verifica qual o vizinho válido
        int vviz = (iy == Ny - 1) ? (i - Nx) : (i + Nx);
        eq.diag = -(2 * isb + 2 * sb + theta);
        eq.vizinhos = {{i - 1, isb}, {i + 1, isb}, {vviz, 2 * sb}};
        eq.rhs = -theta * sistema.bcValues[i];

      } else { // Arestas laterais
        // Verifica qual o vizinho válido
        int hviz = (ix == Nx - 1) ? (i - 1) : (i + 1);
        eq.diag = -(2 * sb + 2 * isb + omega);
        eq.vizinhos = {{i - Nx, sb}, {i + Nx, sb}, {hviz, 2 * isb}};
        eq.rhs = -omega * sistema.bcValues[i];
      }
      return eq;
    }

    // nó interno
    double beta = sistema.config.get_beta();
    eq.diag = -2 * (1 + beta);
    eq.vizinhos = {{i - 1, 1.0}, {i + 1, 1.0}, {i - Nx, beta}, {i + Nx, beta}};
    eq.rhs = 0.0;
    return eq;
  }

public:
  Solver(System &sistema) : sistema(sistema) {}
  // Struct para padronizar o retorno das análises dos solvers
  struct SolverResult {
    double time_ms;
    int iterations;
    double erro_max;
  };

  // Método para montar o sistema linear e resolver por eliminação gaussiana
  SolverResult eliminacao_gaussiana() {
    auto start = std::chrono::high_resolution_clock::now();

    int totalSize = sistema.malha.get_total_size();

    // Monta as matrizes inicializadas com zeros
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(totalSize, totalSize);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(totalSize);

    // Preenche a matriz com os coeficientes calculados pela função
    // monta_equacao()
    for (int i = 0; i < totalSize; i++) {
      Equacao eq = monta_equacao(i);
      A(i, i) = eq.diag;
      for (auto &[j, coef] : eq.vizinhos)
        A(i, j) = coef;
      b(i) = eq.rhs;
    }

    // Eliminação progressiva com pivoteamento
    for (int k = 0; k < totalSize - 1; k++) {

      // Inicia as variáveis com o valor sem pivoteamento
      int max_index = k;
      double max_val = std::abs(A(k, k));

      // Procura o maior valor em módulo na coluna 'k', da linha 'k' para baixo
      for (int i = k + 1; i < totalSize; i++) {
        if (std::abs(A(i, k)) > max_val) {
          max_val = std::abs(A(i, k));
          max_index = i;
        }
      }

      // Se o maior valor encontrado não estiver na linha atual (k), trocamos as
      // linhas
      if (max_index != k) {
        // Troca de linhas da matriz de coeficientes
        A.row(k).swap(A.row(max_index));

        // Troca de linha do vetor de resultados
        std::swap(b(k), b(max_index));
      }

      // Verificação de singularidade/condicionamento
      // Como o determinante é o produtório da diagonal, se um termo for nulo
      // (ou próximo disso) o determinante será nulo e a matriz é singular
      if (std::abs(A(k, k)) < 1e-12) {
        throw std::runtime_error(
            "Sistema singular ou mal condicionado detectado na eliminação!");
      }

      // Eliminação progressiva
      for (int i = k + 1; i < totalSize; i++) {
        double fator = A(i, k) / A(k, k);

        // Zera o elemento abaixo do pivô
        A(i, k) = 0.0;

        // Atualiza o resto da linha
        for (int j = k + 1; j < totalSize; j++) {
          A(i, j) -= fator * A(k, j);
        }
        // Atualiza o vetor de forças 'b'
        b(i) -= fator * b(k);
      }

      // Substituição regressiva para resolver A * T = b
      auto &T = sistema.temperatures; // Referência para atualizar o
                                      // vetor diretamente
      // Loop da última linha da matriz até a primeira
      for (int i = totalSize - 1; i >= 0; --i) {
        double soma = 0.0;
        // Loop que passa de cada coluna j > 1 até a última somando o que já é
        // conhecido
        for (int j = i + 1; j < totalSize; ++j) {
          soma += A(i, j) * T[j];
        }
        // Atualização da temperatura
        T[i] = (b(i) - soma) / A(i, i);
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    return {duration.count(), 1, 0.0};
  }

  SolverResult liebmann(double omega_relax = 1.0) {
    auto start = std::chrono::high_resolution_clock::now();

    int totalSize = sistema.malha.get_total_size();
    std::vector<double> &T = sistema.temperatures; // Referência para atualizar
                                                   // o vetor diretamente

    // Loop por todas as iterações possíveis (para antes em caso de
    // convergência)
    for (int it = 0; it < sistema.config.get_max_It(); it++) {
      double erro_max = 0.0;

      for (int i = 0; i < totalSize; i++) {
        // Da mesma forma que para a matriz, monta a equação para esse nó
        Equacao eq = monta_equacao(i);

        double soma = 0.0;

        // Usando os coeficientes da equação montada calcula a parte
        // relacionada aos vizinhos
        for (auto &[j, coef] : eq.vizinhos)
          soma += coef * T[j];

        // Isola T do nó atual
        double T_novo = (eq.rhs - soma) / eq.diag;

        // Liebmann com relaxamento (default é 1, sem relaxamento)
        T_novo = omega_relax * T_novo + (1 - omega_relax) * T[i];

        // Se o resíduo para o nó atual for maior do que os que já passaram
        // substitui erro_max
        erro_max = std::max(erro_max, std::abs(T_novo - T[i]));

        T[i] = T_novo; // já usa o valor atualizado nos próximos nós
      }

      // Se o sistema converge sai do nó
      if (erro_max < sistema.config.get_tol()) {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;
        return {duration.count(), it + 1, erro_max};
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    return {duration.count(), sistema.config.get_max_It(), 0.0};
  }
};