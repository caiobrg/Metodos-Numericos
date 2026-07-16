// ============================================================================
// Decomposição LU Manual com Substituição Progressiva e Regressiva
// ----------------------------------------------------------------------------
// Resolve um sistema de equações lineares Ax = B (ordem 10) gerado
// aleatoriamente. Implementa a fatoração LU e as substituições manualmente
// (sem resolvedores embutidos da Eigen). Valida os resultados verificando
// se A = LU e se AX ≈ B.
// ============================================================================

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <iostream>

int main() {
  int n = 10; // Ordem do sistema

  // --- Geração aleatória da matriz A (n×n) e do vetor B (n×1) ---
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
  Eigen::VectorXd B = Eigen::VectorXd::Random(n);

  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;
  std::cout << "Matriz A Original:" << std::endl;
  std::cout << A << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;

  // =============================================
  // Decomposição LU Manual (sem pivotamento)
  // =============================================
  // L: triangular inferior com diagonal unitária (L_ii = 1)
  // U: triangular superior (cópia de A, modificada in-place)
  Eigen::MatrixXd L = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd U = A;

  // Eliminação de Gauss para obter L e U
  for (int k = 0; k < n - 1; ++k) {
    // Verifica se o pivô é suficientemente não-nulo
    if (std::abs(U(k, k)) < 1e-12) {
      std::cout
          << "Erro: Divisão por zero detectada na diagonal de U no pivô k = "
          << k << std::endl;
      return 1;
    }
    for (int i = k + 1; i < n; ++i) {
      // Calcula o fator multiplicador e armazena em L
      double fator = U(i, k) / U(k, k);
      L(i, k) = fator;
      // Elimina os elementos abaixo do pivô na coluna k
      for (int j = k; j < n; ++j) {
        U(i, j) -= fator * U(k, j);
      }
    }
  }

  std::cout << "Matriz L (Lower - calculada manualmente):" << std::endl;
  std::cout << L << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;

  std::cout << "Matriz U (Upper - calculada manualmente):" << std::endl;
  std::cout << U << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;

  // =============================================
  // Substituição Progressiva: resolve L * Y = B
  // =============================================
  // Como L tem diagonal unitária, não há divisão por L(i,i)
  Eigen::VectorXd Y(n);
  for (int i = 0; i < n; ++i) {
    double soma = 0.0;
    for (int j = 0; j < i; ++j) {
      soma += L(i, j) * Y(j);
    }
    Y(i) = B(i) - soma; // Y(i) = B(i) - Σ L(i,j)*Y(j) para j < i
  }

  // =============================================
  // Substituição Regressiva: resolve U * X = Y
  // =============================================
  Eigen::VectorXd X(n);
  for (int i = n - 1; i >= 0; --i) {
    double soma = 0.0;
    for (int j = i + 1; j < n; ++j) {
      soma += U(i, j) * X(j);
    }
    // Verifica se o pivô diagonal de U é não-nulo
    if (std::abs(U(i, i)) < 1e-12) {
      std::cout << "Erro: Divisão por zero detectada na diagonal de U na "
                   "substituição regressiva i = "
                << i << std::endl;
      return 1;
    }
    X(i) = (Y(i) - soma) / U(i, i); // X(i) = (Y(i) - Σ U(i,j)*X(j)) / U(i,i)
  }

  std::cout << "Vetor Solução X: " << std::endl << X.transpose() << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;

  // =============================================
  // Prova Real 1: Verifica se A = L * U
  // =============================================
  std::cout << "Prova Real 1 (Verificando se A = L * U):" << std::endl;
  Eigen::MatrixXd LU = L * U;
  std::cout << "Diff (A - LU) norm: " << (A - LU).norm() << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;

  // =============================================
  // Prova Real 2: Verifica se A * X = B
  // =============================================
  Eigen::VectorXd verificando = A * X;
  std::cout << "Prova Real 2 (Verificando se A * X = B):" << std::endl;
  std::cout << "Valores esperados de B: " << B.transpose() << std::endl;
  std::cout << "Valores obtidos:      " << verificando.transpose() << std::endl;
  std::cout << "Diff norm: " << (verificando - B).norm() << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;

  return 0;
}
