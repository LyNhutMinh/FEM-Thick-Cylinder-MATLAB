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
Q4:
<img width="704" height="628" alt="image" src="https://github.com/user-attachments/assets/e826a10a-332e-423f-b923-bbb34db2aa01" />
<img width="701" height="628" alt="image" src="https://github.com/user-attachments/assets/e7bc4adb-0ffa-4b94-be52-57000d2d9e64" />
<img width="701" height="629" alt="image" src="https://github.com/user-attachments/assets/fd8cd071-1de9-4fff-ae45-06767da56424" />
T3:
<img width="703" height="629" alt="image" src="https://github.com/user-attachments/assets/07312989-3112-4f6d-8f6f-e091d66932d8" />
<img width="700" height="629" alt="image" src="https://github.com/user-attachments/assets/10cc66d0-9fa1-4d29-ac17-2997349532cd" />
<img width="705" height="625" alt="image" src="https://github.com/user-attachments/assets/1ed1570f-d57d-4e2f-81a7-413726c4734a" />






## 🎓 Academic Context
This project demonstrates foundational knowledge in computational mechanics, structural analysis, and algorithm translation from mathematical formulations into working code. It is suitable as a reference for mechanical engineering coursework and FEA algorithm studies.
