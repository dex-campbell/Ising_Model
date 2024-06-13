# 2D Ising Model Simulation

This repository contains a Python implementation of the 2D Ising model. The Ising model is a mathematical model used in statistical mechanics to understand phase transitions in ferromagnetic systems.

## Features

- Simulate the 2D Ising model on a lattice
- Calculate physical quantities such as energy, magnetization, specific heat capacity, and susceptibility
- Plot the results of the simulation

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/yourusername/ising-model.git
    ```
2. Navigate to the repository directory:
    ```sh
    cd ising-model
    ```
3. Install the required dependencies:
    ```sh
    pip install -r requirements.txt
    ```

## Usage

1. Run the script:
    ```sh
    python ising3.2.py
    ```
2. The script will prompt you to enter simulation parameters such as the magnetic field strength, temperature range, and lattice size.

3. After entering the parameters, the simulation will run and plot the results.

## Parameters

- **B (float)**: Magnetic field strength
- **T_min (float)**: Minimum temperature for the simulation
- **T_max (float)**: Maximum temperature for the simulation
- **N (int)**: Size of the lattice (NxN)
- **MCS (int)**: Number of Monte Carlo steps for the simulation

## Functions

- `plot_results(Temperature, average_energy, mag, heat_cap, susceptibility_arr, B, N, n)`: Plots the energy, magnetization, specific heat capacity, and susceptibility as functions of temperature.
- `calc_energy(lattice, N, B)`: Calculates the energy of the lattice.
- `calc_mag(lattice, N)`: Calculates the magnetization of the lattice.
- `delta_energy(lattice, i, j, N, B)`: Calculates the change in energy for a given spin configuration.
- `get_parameters()`: Prompts the user to input simulation parameters.

## Example

Here's an example of how to run the script:

```sh
python ising3.2.py
