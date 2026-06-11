#include <cmath>
#include <eigen3/Eigen/Dense>
#include <iostream>

int main() {
  int n = 10;

  // Use Eigen to generate random Matrix A and Vector B
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

  // Manual LU Decomposition using nested for loops
  Eigen::MatrixXd L = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd U = A;

  for (int k = 0; k < n - 1; ++k) {
    if (std::abs(U(k, k)) < 1e-12) {
      std::cout
          << "Erro: Divisão por zero detectada na diagonal de U no pivô k = "
          << k << std::endl;
      return 1;
    }
    for (int i = k + 1; i < n; ++i) {
      double fator = U(i, k) / U(k, k);
      L(i, k) = fator;
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

  // Manual Forward Substitution to solve L * Y = B
  Eigen::VectorXd Y(n);
  for (int i = 0; i < n; ++i) {
    double soma = 0.0;
    for (int j = 0; j < i; ++j) {
      soma += L(i, j) * Y(j);
    }
    Y(i) = B(i) - soma;
  }

  // Manual Backward Substitution to solve U * X = Y
  Eigen::VectorXd X(n);
  for (int i = n - 1; i >= 0; --i) {
    double soma = 0.0;
    for (int j = i + 1; j < n; ++j) {
      soma += U(i, j) * X(j);
    }
    if (std::abs(U(i, i)) < 1e-12) {
      std::cout << "Erro: Divisão por zero detectada na diagonal de U na "
                   "substituição regressiva i = "
                << i << std::endl;
      return 1;
    }
    X(i) = (Y(i) - soma) / U(i, i);
  }

  std::cout << "Vetor Solução X: " << std::endl << X.transpose() << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;

  // Proof 1: L * U = A
  std::cout << "Prova Real 1 (Verificando se A = L * U):" << std::endl;
  Eigen::MatrixXd LU = L * U;
  std::cout << "Diff (A - LU) norm: " << (A - LU).norm() << std::endl;
  std::cout << "---------------------------------------------------------------"
               "-------------------------------------"
            << std::endl;

  // Proof 2: A * X = B
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
