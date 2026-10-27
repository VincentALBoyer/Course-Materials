# Branch-and-Cut Workshop: notebooks

Hands-on notebooks for a 2-hour online workshop on Branch-and-Cut, from the simplex method to the
Traveling Salesperson Problem.

**Instructors:** Prof. Angélica Salazar (Universidad Autónoma de Nuevo León) and
Prof. Vincent Boyer (Tecnológico de Monterrey).
Questions: maria.salazaragl@uanl.edu.mx · vincent.boyer@tec.mx

Two problems run through the workshop: **E**, solved step by step on the slides, and **P**, which you
solve yourself with the same methods. The notebooks do the arithmetic (exact fractions, models written in
PuLP), so you can focus on the concepts.

| Notebook | Topic | Activity |
|----------|-------|----------|
| `01_Simplex_Matrix_Form.ipynb` | The simplex method in matrix form: B, B⁻¹, reduced costs, ratio test, update of B⁻¹ | 2 |
| `02_BranchAndBound_PuLP.ipynb` | Branch-and-Bound: PuLP solves each node, you build and explain the tree | 3 |
| `03_Gomory_Cuts.ipynb` | Gomory cuts from the optimal basis, numerical stability, Chvátal–Gomory cuts | 4 and 5 |
| `04_Knapsack_Cover_Cuts.ipynb` | Knapsack cover inequalities and separation (optional, at home) | – |
| `05_TSP_Branch_and_Cut.ipynb` | TSP: be the separation oracle, then code it | 6 |

## How to run

- **Google Colab:** download a notebook, then *File → Upload notebook* at https://colab.research.google.com
  (or open it from Google Drive with *Open with → Google Colaboratory*).
- **Locally:** Python 3 with Jupyter (or VS Code). The first cell of each notebook installs PuLP 4 with the HiGHS solver if needed;
  the notebooks also use numpy, matplotlib and networkx.

Run the cells from top to bottom. Cells marked **✍️ TODO** are for you.
