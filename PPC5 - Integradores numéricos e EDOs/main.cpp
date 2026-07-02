#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

class ConfiguracaoSimulacao {
private:
  // Atributos encapsulados (privados), não podem ser alterados quando o
  // programa está rodando, definidos na criação do objeto
  double chuteInicial;
  double passoIntegracao;
  double etaMax;
  double tol;
  int itMax;

  // Método auxiliar privado (sobrecarga para double)
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

  // Método auxiliar privado (sobrecarga para int)
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
  // Construtor: Inicializa o objeto com os valores default
  ConfiguracaoSimulacao(double chuteInicial = 0.2,
                        double passoIntegracao = 0.01, double etaMax = 10.0,
                        double tol = 1.0e-5, int itMax = 100)
      : chuteInicial(chuteInicial), passoIntegracao(passoIntegracao),
        etaMax(etaMax), tol(tol), itMax(itMax) {}

  // Método público para iniciar a interação com o usuário e capturar os inputs
  void configurarPeloTerminal() {
    chuteInicial = capturarInput(chuteInicial, "s (chute inicial)");
    passoIntegracao =
        capturarInput(passoIntegracao, "\\delta\\eta (passo de integração)");
    etaMax = capturarInput(etaMax, "\\eta_{max} (ponto máximo da integração)");
    tol = capturarInput(tol, "tol (tolerância)");
    itMax = capturarInput(itMax,
                          "itMax (máximo de iterações para o método do tiro)");
  }

  //  Métodos públicos para acessar os valores de forma segura (somente leitura)
  double getChuteInicial() const { return chuteInicial; }
  double getPassoIntegracao() const { return passoIntegracao; }
  double getEtaMax() const { return etaMax; }
  double getTol() const { return tol; }
  int getItMax() const { return itMax; }
};

// Guarda um estado para cada uma das variáveis y1, y2 e y3 (ex.: derivada,
// coeficientes, etc.)
struct estado {
  double y1;
  double y2;
  double y3;
};

class IntegradorRK4 {
private:
  // Passo de integração (\delta \eta), não pode ser alterado durante a execução
  // do programa
  double h;

  estado calcular_derivadas(const estado &yVet) {
    estado derivadas;
    derivadas.y1 = yVet.y2;
    derivadas.y2 = yVet.y3;
    derivadas.y3 = -0.5 * yVet.y1 * yVet.y3;

    return derivadas;
  }

public:
  // Inicializador da classe
  IntegradorRK4(double passo) : h(passo) {}

  // Metodo que recebe o estado no passo atual e retorna o estado no próximo
  // passo
  estado dar_passo(const estado &yOld) {
    // Os k1's são as derivadas (F_j's) avaliadas no ponto inicial (passo i).
    // k_1^(j) = F_j(eta_i; y_m,i para m = 1,2,3)
    estado k1Vet = calcular_derivadas(yOld);

    // Os k2's são as derivadas avaliadas num ponto intermediário que depende de
    // k1. k_2^(j) = F_j(eta_i + h/2; y_m,i + h/2 * k_1^(m) para m = 1,2,3)
    estado yInterm2 = {.y1 = yOld.y1 + 0.5 * h * k1Vet.y1,
                       .y2 = yOld.y2 + 0.5 * h * k1Vet.y2,
                       .y3 = yOld.y3 + 0.5 * h * k1Vet.y3};
    estado k2Vet = calcular_derivadas(yInterm2);

    // De forma similar aos k2's, os k3's são intermediários dependentes k2.
    // k_3^(j) = F_j(eta_i + h/2; y_m,i + h/2 * k_2^(m) para m = 1,2,3)
    estado yInterm3 = {.y1 = yOld.y1 + 0.5 * h * k2Vet.y1,
                       .y2 = yOld.y2 + 0.5 * h * k2Vet.y2,
                       .y3 = yOld.y3 + 0.5 * h * k2Vet.y3};
    estado k3Vet = calcular_derivadas(yInterm3);

    // Os k4's são avaliados ao "final" do passo, mas ainda dependem de k3
    // k_4^(j) = F_j(eta_i + h; y_m,i + h * k_3^(m), para m = 1,2,3)
    estado yInterm4 = {.y1 = yOld.y1 + h * k3Vet.y1,
                       .y2 = yOld.y2 + h * k3Vet.y2,
                       .y3 = yOld.y3 + h * k3Vet.y3};
    estado k4Vet = calcular_derivadas(yInterm4);

    // Atualização final das variáveis auxiliares y_m
    // y_{m,i+1} = y_{m,i} + (k_{1,i}^(m) + 2 * k_{2,i}^(m) +
    //                                      2 * k_{3,i}^(m) + k_{4,i}^(m))/6
    estado y_new = {
        .y1 = yOld.y1 +
              h * (k1Vet.y1 + 2 * k2Vet.y1 + 2 * k3Vet.y1 + k4Vet.y1) / 6,
        .y2 = yOld.y2 +
              h * (k1Vet.y2 + 2 * k2Vet.y2 + 2 * k3Vet.y2 + k4Vet.y2) / 6,
        .y3 = yOld.y3 +
              h * (k1Vet.y3 + 2 * k2Vet.y3 + 2 * k3Vet.y3 + k4Vet.y3) / 6};
    return y_new;
  }
};

class MetodoDoTiro {
private:
  // Parâmetros iniciais para a integração
  ConfiguracaoSimulacao config;

  // Instância do integrador a ser usada pelo método do tiro
  IntegradorRK4 integrador{config.getPassoIntegracao()};

  // Número de pontos que o integrador vai percorrer
  std::size_t N = static_cast<std::size_t>(
      std::round(config.getEtaMax() / config.getPassoIntegracao()) + 1);

  // Vetores que guardam os valores de y ao longo dos 2 chutes iniciais
  std::vector<estado> yVec1;
  std::vector<estado> yVec2;

public:
  // Inicialização da classe com os argumentos de entrada
  MetodoDoTiro(ConfiguracaoSimulacao configuracao) : config(configuracao) {}

  // // Helper function para salvar os resultados de um vetor num arquivo .csv
  void _save_results(const std::vector<estado> &vec, const char *filename,
                     int precision = 15) {

    // Abre o arquivo
    std::string filenameStr = filename;
    std::ofstream outFile("./results/" + filenameStr + ".csv");

    // Verifica se foi possível abrir o arquivo ou para o programa
    if (!outFile.is_open())
      throw std::runtime_error("Não foi possível abrir o arquivo " +
                               filenameStr + " para gravação");

    // Define o número de casas decimais no arquivo
    outFile << std::fixed << std::setprecision(precision);

    // Printa o cabeçalho do arquivo
    outFile << "eta,y1,y2,y3" << std::endl;

    // Printa os dados do vetor
    for (size_t i = 0; i < vec.size(); i++) {
      outFile << i * config.getPassoIntegracao() << "," << vec[i].y1 << ","
              << vec[i].y2 << "," << vec[i].y3 << std::endl;
    }

    // Fecha o arquivo
    outFile.close();
  }

  // Helper function para printar uma linha de iteração para o terminal
  void _print_to_terminal(const std::vector<estado> &vec, int iter) {
    std::cout << std::setw(15) << iter << std::setw(15) << vec[0].y3
              << std::setw(15) << vec[N - 1].y2 << std::setw(15)
              << std::abs(vec[N - 1].y2 - 1) << std::endl;
  }

  // Método principal de execução, executa o método do tiro e salva, pelo
  // menos, o resultado final
  void execucao(bool verbose = false, bool save_all = false) {
    // Variaveis de convergencia
    double err1 = 1.0;
    double err2 = 1.0;
    double tol = config.getTol();
    int iter = 0;

    // Cria o diretório de resultados caso ele não exista
    std::filesystem::create_directory("./results");

    if (verbose) {
      std::cout << std::fixed << std::setprecision(6);
      std::cout << "Iniciando o Metodo do Tiro...\n";
    }

    // Inicializando os vetores para a primeira run
    yVec1.assign(N, {.y1 = 0.0, .y2 = 0.0, .y3 = config.getChuteInicial()});
    yVec2.assign(N, {.y1 = 0.0, .y2 = 0.0, .y3 = config.getChuteInicial() * 2});

    // Passa por cada ponto da malha integrando pelo método RK4
    for (size_t i = 0; i < N - 1; i++) {
      yVec1[i + 1] = integrador.dar_passo(yVec1[i]);
      yVec2[i + 1] = integrador.dar_passo(yVec2[i]);
    }

    // Salva os dados no arquivo de resultados caso isso tenha sido solicitado
    if (save_all) {
      _save_results(yVec1, "iter1");
      _save_results(yVec2, "iter2");
    }

    // Atualização do número de iterações
    iter += 2;

    // Loga os dados no terminal caso o usuário tenha selecionado essa opção
    if (verbose) {
      std::cout << std::setw(15) << "Iter" << std::setw(15) << "Chute"
                << std::setw(15) << "y2(eta_max)" << std::setw(15) << "err"
                << std::endl;
      _print_to_terminal(yVec1, 1);
      _print_to_terminal(yVec2, 2);
    }

    // O erro é a diferença entre o valor de f' em eta_max e a condição de
    // contorno no infinito (f' = 1). Erro para o primeiro chute:
    err1 = yVec1[N - 1].y2 - 1;
    // Erro para o segundo chute
    err2 = yVec2[N - 1].y2 - 1;

    while (std::max(std::abs(err2), std::abs(err1)) >= tol &&
           iter <= config.getItMax()) {

      // Calcula pelo método da secante um novo valor tentativo inicial
      double valor_inicial_new =
          yVec2[0].y3 - err2 * (yVec2[0].y3 - yVec1[0].y3) / (err2 - err1);

      // Na próxima iteração será utilizado o ponto 2 (agora ponto 1) e o ponto
      // 3 (futuramente ponto 2)
      yVec1 = yVec2;
      err1 = err2;
      yVec2.assign(N, {.y1 = 0, .y2 = 0, .y3 = valor_inicial_new});

      // Integra a curva usando o novo chute
      for (size_t i = 0; i < N - 1; i++) {
        yVec2[i + 1] = integrador.dar_passo(yVec2[i]);
      }

      // Erro para o novo valor inicial
      err2 = yVec2[N - 1].y2 - 1;

      iter++;

      if (verbose)
        _print_to_terminal(yVec2, iter);

      if (save_all)
        _save_results(yVec2, ("iter" + std::to_string(iter)).c_str());
    }

    _save_results(yVec2, ("iter" + std::to_string(iter)).c_str());

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "Valor de f\" para convergência: f\"(0) = "
              << std::setprecision(15) << yVec2[0].y3 << std::endl;
    std::cout << "Erro final obtido: |f'(eta_max) - 1| = "
              << std::setprecision(15) << std::abs(yVec2[N - 1].y2 - 1)
              << std::endl;
    std::cout << "Valor obtido para a velocidade adimensional: f'(eta_max) = "
              << std::setprecision(15) << yVec2[N - 1].y2 << std::endl;
    std::cout << "Número de iterações do método do tiro: iter = " << iter
              << std::endl;

    // ========================================================================
    // ANÁLISES FÍSICAS DA CAMADA LIMITE
    // ========================================================================

    std::cout << std::string(80, '-') << std::endl;
    std::cout << "ANALISES DA CAMADA LIMITE\n";

    // 1. Coeficiente de Atrito Local
    // A fórmula é Cf = (2 * f"(0)) / sqrt(Re_x).
    // Isolando a constante numérica, temos Cf * sqrt(Re_x) = 2 * f"(0)
    double constante_Cf = 2.0 * yVec2[0].y3;
    std::cout << "Constante do Coeficiente de Atrito [ Cf * sqrt(Re_x) ] = "
              << constante_Cf << std::endl;

    // 2. Cálculo da espessura da Camada Limite (eta_99) via Interpolação Linear
    double eta99 = 0.0;
    double h = config.getPassoIntegracao();

    for (size_t i = 1; i < N; i++) {
      if (yVec2[i].y2 >= 0.99) {
        double eta_ant = (i - 1) * h;
        double eta_atual = i * h;
        double f_linha_ant = yVec2[i - 1].y2;
        double f_linha_atual = yVec2[i].y2;

        // Fórmula da interpolação linear: y = y0 + (x - x0) * (y1 - y0) / (x1 -
        // x0) Resolvendo para x=0.99:
        eta99 = eta_ant + (0.99 - f_linha_ant) * (eta_atual - eta_ant) /
                              (f_linha_atual - f_linha_ant);

        break; // Encontrou o ponto, para o loop
      }
    }

    std::cout << "Posicao adimensional (eta_99) = " << eta99 << std::endl;
    std::cout << "Coeficiente da espessura (C_delta) = " << eta99 << std::endl;
    std::cout << std::string(80, '-') << std::endl;
  }
};

int main() {
  // Instancia o objeto de configuração
  ConfiguracaoSimulacao config;

  // Pede para o usuário preencher/confirmar os dados, pode ser desativado
  // comentando a linha, valor default serão usados
  config.configurarPeloTerminal();

  std::cout << "\n--- Resumo da Configuracao ---\n";
  std::cout << "Chute Inicial: " << config.getChuteInicial() << "\n";
  std::cout << "Passo: " << config.getPassoIntegracao() << "\n";
  std::cout << "Iteracoes Maximas: " << config.getItMax() << "\n";

  // Cria uma instância da classe MetodoDoTiro
  MetodoDoTiro metodo{config};

  // Inicia o loop de execucao principal
  // Opções .execucao(verbose = false, save_all = false)
  metodo.execucao(true, true);
  return 0;
}
