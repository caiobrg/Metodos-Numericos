#include <bits/stdc++.h>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// Vetor de estados
struct estado {
  double x;
  double y;
};

class configuracao_simulacao {
private:
  // Variáveis privadas não podem ser alteradas durante a execução do programa

  // Constantes físicas do problemas
  double G;
  double m1;
  double m2;

  // Valores iniciais para os vetores r
  estado r1_ini;
  estado r2_ini;

  // Valores iniciais da velocidade
  estado v1_ini;
  estado v2_ini;

  double passoIntegracao;
  double tempoTotal;

  std::string metodo;

  // Método auxiliar privado (sobrecarga para string)
  std::string capturarInput(std::string valorPadrao, const std::string &nome) {
    std::string input;
    std::cout << "Digite o valor de " << nome << " [default: " << valorPadrao
              << "]: ";
    std::getline(std::cin, input);

    // Transforma o input para lowercase
    std::transform(input.begin(), input.end(), input.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (!input.empty()) {
      return input;
    }
    return valorPadrao;
  }

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
  configuracao_simulacao(std::string metodo = "euler",
                         double passoIntegracao = 1.0e-3,
                         double tempoTotal = 10, double G = 1, double m1 = 1,
                         double m2 = 1, estado r1_ini = {.x = -0.5, .y = 0},
                         estado r2_ini = {.x = 0.5, .y = 0},
                         estado v1_ini = {.x = 0, .y = -0.5},
                         estado v2_ini = {.x = 0, .y = 0.5})
      : passoIntegracao(passoIntegracao), tempoTotal(tempoTotal), G(G), m1(m1),
        m2(m2), r1_ini(r1_ini), r2_ini(r2_ini), v1_ini(v1_ini), v2_ini(v2_ini),
        metodo(metodo) {}

  void configurar_pelo_terminal() {
    passoIntegracao =
        capturarInput(passoIntegracao, "\\delta t (passo de integração)");
    tempoTotal = capturarInput(tempoTotal, "t_max (Tempo total de simulação)");
    G = capturarInput(G, "G (Constante gravitacional)");
    m1 = capturarInput(m1, "m1 (massa do corpo 1)");
    m2 = capturarInput(m2, "m2 (massa do corpo 2)");
    r1_ini.x =
        capturarInput(r1_ini.x, "r1_x(0) (posição inicial x do corpo 1)");
    r1_ini.y =
        capturarInput(r1_ini.y, "r1_y(0) (posição inicial y do corpo 1)");
    r2_ini.x =
        capturarInput(r2_ini.x, "r2_x(0) (posição inicial x do corpo 2)");
    r2_ini.y =
        capturarInput(r2_ini.y, "r2_y(0) (posição inicial y do corpo 2)");
    v1_ini.x =
        capturarInput(v1_ini.x, "v1_x(0) (velocidade inicial x do corpo 1)");
    v1_ini.y =
        capturarInput(v1_ini.y, "v1_y(0) (velocidade inicial y do corpo 1)");
    v2_ini.x =
        capturarInput(v2_ini.x, "v2_x(0) (velocidade inicial x do corpo 2)");
    v2_ini.y =
        capturarInput(v2_ini.y, "v2_y(0) (velocidade inicial y do corpo 2)");
    while (true) {
      metodo = capturarInput("euler", "metodo (escolha 'euler' ou 'leapfrog')");
      if (metodo == "euler" || metodo == "leapfrog") {
        break;
      }
      std::cout << "Método inválido! Por favor, escolha 'euler' ou 'leapfrog'."
                << std::endl;
    }
  }

  double get_passo_integracao() const { return passoIntegracao; }
  double get_tempo_total() const { return tempoTotal; }
  double get_G() const { return G; }
  double get_m1() const { return m1; }
  double get_m2() const { return m2; }
  estado get_r1_ini() const { return r1_ini; }
  estado get_r2_ini() const { return r2_ini; }
  estado get_v1_ini() const { return v1_ini; }
  estado get_v2_ini() const { return v2_ini; }
  std::string get_metodo() const { return metodo; }
};

class integrador {
private:
  configuracao_simulacao config;

  std::string metodo = config.get_metodo();

  double dt = config.get_passo_integracao();
  double G = config.get_G();
  double m1 = config.get_m1();
  double m2 = config.get_m2();

  std::vector<estado> calcular_derivadas_v(const estado &r1, const estado &r2) {
    estado v1;
    estado v2;

    double alpha =
        std::sqrt(std::pow((r2.x - r1.x), 2) + std::pow((r2.y - r1.y), 2));

    v1.x = G * m2 * (r2.x - r1.x) / std::pow(alpha, 3);
    v1.y = G * m2 * (r2.y - r1.y) / std::pow(alpha, 3);
    v2.x = -G * m1 * (r2.x - r1.x) / std::pow(alpha, 3);
    v2.y = -G * m1 * (r2.y - r1.y) / std::pow(alpha, 3);

    std::vector<estado> result(2);

    result[0] = v1;
    result[1] = v2;

    return result;
  }

public:
  integrador(configuracao_simulacao config) : config(config) {}

  std::vector<estado> dar_passo(const estado &r1_old, const estado &v1_old,
                                const estado &r2_old, const estado &v2_old) {

    estado r1_new;
    estado r2_new;
    estado v1_new;
    estado v2_new;

    std::vector<estado> acc = calcular_derivadas_v(r1_old, r2_old);

    if (metodo == "euler") {
      // Atualiza as velocidades usando a aceleração calculada a partir da
      // posição no passo anterior
      v1_new.x = v1_old.x + dt * acc[0].x;
      v1_new.y = v1_old.y + dt * acc[0].y;
      v2_new.x = v2_old.x + dt * acc[1].x;
      v2_new.y = v2_old.y + dt * acc[1].y;

      // Atualiza a posição usando a velocidade calculada no passo anterior
      r1_new.x = r1_old.x + dt * v1_old.x;
      r1_new.y = r1_old.y + dt * v1_old.y;
      r2_new.x = r2_old.x + dt * v2_old.x;
      r2_new.y = r2_old.y + dt * v2_old.y;

    } else if (metodo == "leapfrog") {

      estado v1_interm;
      estado v2_interm;

      // Calcula a velocidade num passo intermediário
      v1_interm.x = v1_old.x + (dt / 2) * acc[0].x;
      v1_interm.y = v1_old.y + (dt / 2) * acc[0].y;
      v2_interm.x = v2_old.x + (dt / 2) * acc[1].x;
      v2_interm.y = v2_old.y + (dt / 2) * acc[1].y;

      // Atualiza os raios usando a velocidade nesse passo intermediário
      r1_new.x = r1_old.x + dt * v1_interm.x;
      r1_new.y = r1_old.y + dt * v1_interm.y;
      r2_new.x = r2_old.x + dt * v2_interm.x;
      r2_new.y = r2_old.y + dt * v2_interm.y;

      // Atualiza a aceleração usando as novas posições
      acc = calcular_derivadas_v(r1_new, r2_new);

      // Por fim, atualiza as velocidades usando o ponto médio e a nova
      // aceleração
      v1_new.x = v1_interm.x + (dt / 2) * acc[0].x;
      v1_new.y = v1_interm.y + (dt / 2) * acc[0].y;
      v2_new.x = v2_interm.x + (dt / 2) * acc[1].x;
      v2_new.y = v2_interm.y + (dt / 2) * acc[1].y;
    }

    std::vector<estado> results(4);
    // Salva os novos resultados num vetor de resultados
    results[0] = r1_new;
    results[1] = v1_new;
    results[2] = r2_new;
    results[3] = v2_new;

    return results;
  }
};

class solver {
private:
  configuracao_simulacao config;
  integrador integrad{config};

  std::size_t N = static_cast<std::size_t>(
      std::round(config.get_tempo_total() / config.get_passo_integracao()) + 1);

  // Vetores que guardam os valores ao longo da integração
  std::vector<estado> r1Vec;
  std::vector<estado> v1Vec;
  std::vector<estado> r2Vec;
  std::vector<estado> v2Vec;

public:
  solver(configuracao_simulacao config) : config(config) {}

  void _save_results(const std::vector<estado> &r1Vec,
                     const std::vector<estado> &v1Vec,
                     const std::vector<estado> &r2Vec,
                     const std::vector<estado> &v2Vec,
                     const std::string &filename, int precision = 15) {
    // Abre o arquivo
    std::ofstream outFile("./results/" + filename + ".csv");

    // Verifica se foi possível abrir o arquivo ou para o programa
    if (!outFile.is_open())
      throw std::runtime_error("Não foi possível abrir o arquivo " + filename +
                               " para gravação");
    // Define o número de casas decimais no arquivo
    outFile << std::fixed << std::setprecision(precision);

    // Printa o cabeçalho do arquivo
    outFile << "t,x1,y1,vx1,vy1,x2,y2,vx2,vy2,E" << std::endl;

    // Printa os dados do vetor
    for (size_t i = 0; i < r1Vec.size(); i++) {
      double xRes = r2Vec[i].x - r1Vec[i].x;
      double yRes = r2Vec[i].y - r1Vec[i].y;

      double modRes = std::sqrt(std::pow(xRes, 2) + std::pow(yRes, 2));

      double modv1 =
          std::sqrt(v1Vec[i].x * v1Vec[i].x + v1Vec[i].y * v1Vec[i].y);
      double modv2 =
          std::sqrt(v2Vec[i].x * v2Vec[i].x + v2Vec[i].y * v2Vec[i].y);

      double K = 0.5 * config.get_m1() * std::pow(modv1, 2) +
                 0.5 * config.get_m2() * std::pow(modv2, 2);

      double U = -config.get_G() * (config.get_m1() * config.get_m2()) / modRes;

      double E = K + U;

      outFile << i * config.get_passo_integracao() << "," << r1Vec[i].x << ","
              << r1Vec[i].y << "," << v1Vec[i].x << "," << v1Vec[i].y << ","
              << r2Vec[i].x << "," << r2Vec[i].y << "," << v2Vec[i].x << ","
              << v2Vec[i].y << "," << E << std::endl;
    }

    // Fecha o arquivo
    outFile.close();
  }

  void execucao() {

    // Cria o diretório de resultados caso ele não exista
    std::filesystem::create_directory("./results");

    r1Vec.assign(N, config.get_r1_ini());
    r2Vec.assign(N, config.get_r2_ini());
    v1Vec.assign(N, config.get_v1_ini());
    v2Vec.assign(N, config.get_v2_ini());

    for (size_t i = 0; i < N - 1; i++) {
      // Vetor temporário para armazenar os resultados de uma etapa de
      // integração
      std::vector<estado> vetTemp(4);

      // Realiza uma etapa de integração
      vetTemp = integrad.dar_passo(r1Vec[i], v1Vec[i], r2Vec[i], v2Vec[i]);

      // Atualiza os valores dos vetores
      r1Vec[i + 1] = vetTemp[0];
      v1Vec[i + 1] = vetTemp[1];
      r2Vec[i + 1] = vetTemp[2];
      v2Vec[i + 1] = vetTemp[3];
    }

    _save_results(r1Vec, v1Vec, r2Vec, v2Vec,
                  "resultados_" + config.get_metodo());
  }
};

int main() {
  configuracao_simulacao configEuler{"euler"};
  solver solverEuler{configEuler};

  solverEuler.execucao();

  configuracao_simulacao configLeapfrog{"leapfrog"};
  solver solver{configLeapfrog};

  solver.execucao();

  return 1;
}
