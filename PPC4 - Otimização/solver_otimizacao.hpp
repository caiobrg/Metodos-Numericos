#pragma once
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

/**
 * @brief Função objetivo a ser maximizada: f(x,y) = 2xy + 2x - x^2 - 2y^2
 */
inline double f(double x, double y) {
  return 2.0 * x * y + 2.0 * x - x * x - 2.0 * y * y;
}

/**
 * @brief Componente x do gradiente da função objetivo.
 */
inline double df_dx(double x, double y) { return 2.0 * y + 2.0 - 2.0 * x; }

/**
 * @brief Componente y do gradiente da função objetivo.
 */
inline double df_dy(double x, double y) { return 2.0 * x - 4.0 * y; }

/**
 * @brief Calcula o módulo do vetor gradiente (erro).
 */
inline double grad_norm(double dx, double dy) {
  return std::sqrt(dx * dx + dy * dy);
}

/**
 * @brief Função paramétrica g(h) ao longo da direção de busca.
 */
inline double g_param(double x, double y, double px, double py, double h) {
  return f(x + h * px, y + h * py);
}

/**
 * @brief Realiza a busca unidimensional por interpolação quadrática de três
 * pontos.
 * @param x Posição atual x.
 * @param y Posição atual y.
 * @param px Direção de busca em x.
 * @param py Direção de busca em y.
 * @return Passo ótimo h*.
 */
inline double line_search(double x, double y, double px, double py) {
  // Escolha inicial de três pontos para avaliação do passo
  double h0 = 0.0;
  double h1 = 0.1;
  double h2 = 0.2;

  double g0 = g_param(x, y, px, py, h0);
  double g1 = g_param(x, y, px, py, h1);
  double g2 = g_param(x, y, px, py, h2);

  double num = (h0 * h0 - h2 * h2) * g1 + (h2 * h2 - h1 * h1) * g0 +
               (h1 * h1 - h0 * h0) * g2;
  double den = (h0 - h2) * g1 + (h2 - h1) * g0 + (h1 - h0) * g2;

  double h_star = 0.0;

  // Fallback: se o denominador for nulo ou muito pequeno
  if (std::abs(den) < 1e-12) {
    if (g0 >= g1 && g0 >= g2)
      h_star = h0;
    else if (g1 >= g0 && g1 >= g2)
      h_star = h1;
    else
      h_star = h2;
  } else {
    h_star = 0.5 * (num / den);

    // Garantir que a interpolação forneceu um ponto que de fato maximiza a
    // função, caso contrário, utiliza-se o melhor
    // ponto inicial.
    double g_star = g_param(x, y, px, py, h_star);
    if (g0 > g_star || g1 > g_star || g2 > g_star) {
      if (g0 >= g1 && g0 >= g2)
        h_star = h0;
      else if (g1 >= g0 && g1 >= g2)
        h_star = h1;
      else
        h_star = h2;
    }
  }

  return h_star;
}

/**
 * @brief Implementa o método de otimização do Aclive Máximo (Steepest Ascent).
 */
inline void steepest_ascent(double x0, double y0, double tol, int max_iter,
                            const std::string &filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Erro ao abrir o arquivo: " << filename << "\n";
    return;
  }

  file << std::fixed << std::setprecision(6);

  double x = x0;
  double y = y0;
  int iter = 0;
  double h = 0.0;

  double dx = df_dx(x, y);
  double dy = df_dy(x, y);
  double erro = grad_norm(dx, dy);

  std::cout << "\n--- Metodo do Aclive Maximo ---\n";

  while (erro > tol && iter < max_iter) {
    std::cout << "Iter: " << iter << " | x: " << x << " | y: " << y
              << " | f: " << f(x, y) << " | h: " << h << " | erro: " << erro
              << "\n";

    file << iter << " " << erro << " " << h << " " << x << " " << y << " " << dx
         << " " << dy << "\n";

    // Direção de busca para maximização é o próprio gradiente
    double px = dx;
    double py = dy;

    // Otimização unidimensional para obtenção do passo
    h = line_search(x, y, px, py);

    // Atualização da posição
    x += h * px;
    y += h * py;
    iter++;

    // Atualização do gradiente e cálculo do novo erro
    dx = df_dx(x, y);
    dy = df_dy(x, y);
    erro = grad_norm(dx, dy);
  }

  // Imprime o ponto de convergência
  std::cout << "Iter: " << iter << " | x: " << x << " | y: " << y
            << " | f: " << f(x, y) << " | h: " << h << " | erro: " << erro
            << "\n";
  file << iter << " " << erro << " " << h << " " << x << " " << y << " " << dx
       << " " << dy << "\n";

  file.close();
}

/**
 * @brief Implementa o método de Gradientes Conjugados (Fletcher-Reeves).
 */
inline void conjugate_gradients(double x0, double y0, double tol, int max_iter,
                                const std::string &filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Erro ao abrir o arquivo: " << filename << "\n";
    return;
  }

  file << std::fixed << std::setprecision(6);

  double x = x0;
  double y = y0;
  int iter = 0;
  double h = 0.0;

  double dx = df_dx(x, y);
  double dy = df_dy(x, y);
  double erro = grad_norm(dx, dy);

  // Na primeira iteração, a direção de busca é o próprio gradiente
  double px = dx;
  double py = dy;

  std::cout << "\n--- Metodo dos Gradientes Conjugados (Fletcher-Reeves) ---\n";

  while (erro > tol && iter < max_iter) {
    std::cout << "Iter: " << iter << " | x: " << x << " | y: " << y
              << " | f: " << f(x, y) << " | h: " << h << " | erro: " << erro
              << "\n";

    file << iter << " " << erro << " " << h << " " << x << " " << y << " " << dx
         << " " << dy << "\n";

    h = line_search(x, y, px, py);

    x += h * px;
    y += h * py;

    double dx_new = df_dx(x, y);
    double dy_new = df_dy(x, y);
    double erro_new = grad_norm(dx_new, dy_new);

    if (erro_new <= tol) {
      iter++;
      dx = dx_new;
      dy = dy_new;
      erro = erro_new;
      break;
    }

    // Atualização da direção via parâmetro de Fletcher-Reeves
    double beta = (dx_new * dx_new + dy_new * dy_new) / (dx * dx + dy * dy);

    px = dx_new + beta * px;
    py = dy_new + beta * py;

    dx = dx_new;
    dy = dy_new;
    erro = erro_new;

    iter++;
  }

  // Imprime o ponto de convergência
  std::cout << "Iter: " << iter << " | x: " << x << " | y: " << y
            << " | f: " << f(x, y) << " | h: " << h << " | erro: " << erro
            << "\n";
  file << iter << " " << erro << " " << h << " " << x << " " << y << " " << dx
       << " " << dy << "\n";

  file.close();
}