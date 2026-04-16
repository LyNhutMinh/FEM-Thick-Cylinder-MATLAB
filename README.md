# FEM-Thick-Cylinder-MATLAB
Eng: A MATLAB-based Finite Element Analysis (FEA) solver for thick-walled cylinders under internal/external pressure (Plane Strain). Features Q4/T3 mesh generation, mesh convergence study, and comparison with exact Lame's analytical solution.

Vie: Chương trình MATLAB giải bài toán ống trụ dày chịu áp suất (Plane Strain) bằng phương pháp phần tử hữu hạn (FEA). Hỗ trợ phần tử Q4/T3, đánh giá hội tụ lưới và so sánh với nghiệm giải tích Lame.

# Finite Element Analysis of Thick-Walled Cylinder (Plane Strain)

## 📌 Project Overview
This repository contains a purely MATLAB-based Finite Element Analysis (FEA) solver developed to analyze thick-walled cylinders subjected to internal and external pressures. The problem is modeled under the **Plane Strain** assumption. 

The program not only calculates numerical solutions using custom-built stiffness matrices but also automatically compares the FEM results against **Lame's exact analytical solutions** to evaluate accuracy.

## ✨ Key Features
* **Custom Mesh Generation:** Procedurally generates structured meshes for $1/4$ symmetry models and full 360-degree models.
* **Element Types:** Supports both 4-node Quadrilateral (Q4) and 3-node Constant Strain Triangle (CST / T3) elements.
* **Mesh Convergence Study:** Automatically iterates through different mesh densities to plot displacement error convergence.
* **Comprehensive Post-Processing:** * Calculates and visualizes radial displacement ($u_r$).
  * Computes radial stress ($\sigma_r$) and hoop stress ($\sigma_\theta$).
  * Evaluates equivalent **Von Mises stress** for yield criterion checking.
* **Analytical Comparison:** Plots direct comparisons between FEM numerical outputs and Lame's exact equations.

## 🛠️ Tech Stack & Methodology
* **Language:** MATLAB
* **Formulation:** 2D Plane Strain (Biến dạng phẳng)
* **Integration Method:** $2 \times 2$ Gauss Quadrature for Q4 elements.
* **Matrix Operations:** Utilizes sparse matrix assembly (Triplet format) for highly optimized and fast computational performance.

## 🚀 How to Run
1. Clone the repository to your local machine.
2. Open MATLAB and run the main script `BTLFEMcubecode.m`.
3. An input dialog will prompt you to enter material properties and boundary conditions:
   * Young's Modulus ($E$) & Poisson's ratio ($\nu$)
   * Inner/Outer Radius ($R_i$, $R_o$)
   * Internal/External Pressure ($p_i$, $p_o$)
   * Desired mesh sizing and element type (Q4 or T3)
4. The script will output the convergence plot, model mesh, and contour plots of stresses and displacements.

## 📊 Sample Outputs


## 🎓 Academic Context
This project demonstrates foundational knowledge in computational mechanics, structural analysis, and algorithm translation from mathematical formulations into working code. It is suitable as a reference for mechanical engineering coursework and FEA algorithm studies.
