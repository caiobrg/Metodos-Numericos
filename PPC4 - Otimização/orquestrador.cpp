#include "solver_otimizacao.hpp"
#include <iostream>
#include <ostream>

int main() {
  std::cout << "========================================\n";
  std::cout << "  Orquestrador - Metodos de Otimizacao  \n";
  std::cout << "========================================\n\n";

  double x0, y0;
  std::cout << "Insira o valor inicial de x0: ";
  std::cin >> x0;
  std::cout << "Insira o valor inicial de y0: ";
  std::cin >> y0;

  double tol = 1e-6;
  int max_iter = 1000;

  const char *filename1 = "output1.dat";
  const char *filename2 = "output2.dat";

  // Resolve utilizando a estratégia de Aclive Máximo
  steepest_ascent(x0, y0, tol, max_iter, filename1);

  // Resolve utilizando a estratégia de Gradientes Conjugados (Fletcher-Reeves)
  conjugate_gradients(x0, y0, tol, max_iter, filename2);

  std::cout << "\n========================================\n";
  std::cout << "Simulacoes concluidas com sucesso!\n";
  std::cout << "-> Aclive Maximo salvo em: " << filename1 << std::endl;
  std::cout << "-> Gradientes Conjugados salvo em: " << filename2 << std::endl;
  std::cout << "-> Execute o script Python para visualizar as trajetorias e "
               "curvas de nivel.\n";
  std::cout << "========================================\n";

  return 0;
}