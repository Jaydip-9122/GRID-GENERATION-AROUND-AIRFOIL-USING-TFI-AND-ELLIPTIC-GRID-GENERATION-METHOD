# Airfoil Grid Generation Toolkit (SOLVER.py)

This script generates computational grids around an airfoil geometry using **Transfinite Interpolation (TFI)** and **Elliptic Grid Generation** methods.  
It also supports spline smoothing of the airfoil surface coordinates and visualizes the generated grids.

---

## Features

- **Read Airfoil Data**: Reads raw airfoil coordinates from CSV (`NACA63412coordinates.csv`).
- **Spline Smoothing**: Uses cubic spline interpolation for smooth airfoil surfaces.
- **TFI Grid Generation**: Quickly generates a structured grid based on boundary interpolation.
- **Elliptic Grid Generation**: Refines the grid using iterative elliptic PDE solver for smoothness.
- **Visualization**: Plots grids and airfoil using Matplotlib.

---

## Requirements

- Python 3.x  
- Required Libraries:  
  - `numpy`  
  - `pandas`  
  - `matplotlib`  
  - `scipy`

Install dependencies:
```bash
pip install numpy pandas matplotlib scipy
```

---

## Files Required

- **SOLVER.py** (this script)  
- **NACA63412coordinates.csv** (airfoil coordinates file with header and 2D points)

---

## How to Run

1. Place the script and `NACA63412coordinates.csv` in the same folder.
2. Open terminal in the folder and run:
```bash
python SOLVER.py
```
3. The script will:
   - Read airfoil data
   - Generate a grid using **TFI**
   - Generate a refined grid using **Elliptic Solver**
   - Display both grids with Matplotlib plots

---

## Code Structure

### Key Functions:
- `NACA_AIRFOIL(file_path)` → Reads and formats airfoil coordinates.
- `CUBICSPLINE_FITTING(x, y, n_points)` → Smooths coordinates.
- `TFI_GRID_GENERATION(NX, NY, chord_length, radius)` → Generates initial TFI grid.
- `ELLIPTIC_GRID_GENERATION(X, Y, iterations, tol)` → Improves grid smoothness.
- `plot_grid(X, Y, x_airfoil, y_airfoil, title)` → Displays grid plot.

### Main Execution:
- Grid size: `NX = 60`, `NY = 60`
- Chord length = `1`
- Outer boundary radius = `10`
- Displays:
  - **TFI Grid**
  - **Elliptic Grid**

---

## Example Output (CLI)
```
Iteration 0, Error: 0.05432
Iteration 1, Error: 0.02641
...
Iteration 49, Error: 0.00092
```
Plots will pop up showing the generated structured grids.

---

## Notes
- The airfoil CSV file must contain coordinate points with at least 5 header rows skipped (`skiprows=5` in code).
- Output grids are suitable for CFD preprocessing or academic visualization.
- Parameters like `NX`, `NY`, and outer radius can be modified inside the script.

