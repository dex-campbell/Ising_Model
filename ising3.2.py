import matplotlib.pyplot as plt
import numpy as np
import random
import time
from tqdm import tqdm

# Physical constants
boltz_cont = 1
bohr_mag = 1


# Function to plot the results
def plot_results(Temperature,average_energy,mag,heat_cap,susceptibility_arr, B, N, n):
    f = plt.figure(figsize=(18, 10))

    plots = f.add_subplot(2, 2, 1)
    plt.plot(Temperature, average_energy, marker='.', color='Red')
    plt.xlabel("Temperature (T)", fontsize=10)
    plt.ylabel("Energy ", fontsize=10)
    plt.axis('tight')
    
    plots = f.add_subplot(2, 2, 2)
    plt.plot(Temperature, heat_cap, marker='.', color='Black')
    plt.xlabel("Temperature (T)", fontsize=10)
    plt.ylabel("Specific Heat Capacity", fontsize=10)
    plt.axis('tight')
    
    plots = f.add_subplot(2, 2, 3)
    plt.plot(Temperature, abs(mag), marker='.', color='Blue')
    plt.xlabel("Temperature (T)", fontsize=10)
    plt.ylabel("Magnetization ", fontsize=10)
    plt.axis('tight')

    plots = f.add_subplot(2, 2, 4)
    plt.plot(Temperature, susceptibility_arr, marker='.', color='Green')
    plt.xlabel("Temperature (T)", fontsize=10)
    plt.ylabel("Susceptibility ", fontsize=10)
    plt.axis('tight')

    plt.show()

# Energy of the lattice calculations
def calc_energy(lattice, N, B):
    energy = 0
    for i in range(N):
        for j in range(N):
            energy += -delta_energy(lattice, i, j, N, B)/2
    return energy/N**2

# Calculate the Magnetization of a given configuration
def calc_mag(lattice, N):
    mag = np.abs(lattice.sum())
    return mag/N**2

# Calculate interaction energy between spins
def delta_energy(lattice, i, j, N, B):
    # Periodic boundaries
    left = lattice[i,(j-1)%N]
    right = lattice[i,(j+1)%N]
    top = lattice[(i-1)%N,j]
    bottom = lattice[(i+1)%N,j]

    energy = 2*lattice[i, j]*(right+left+top+bottom ) - bohr_mag*B*lattice.sum()
    return energy
# Function to get user input for simulation parameters
def get_parameters():
    B = float(input("Enter the magnetic field strength (B): "))  # Prompt user for magnetic field strength
    N = int(input("Enter the lattice size (width): "))  # Prompt user for lattice size
    n = int(input("Enter the number of times to run the algorithm: "))  # Prompt user for number of Monte Carlo steps
    single_temp = input("Do you want to simulate at a single temperature? (yes/no): ").lower().strip()  # Prompt user for single or multiple temperatures
    if single_temp == 'yes':
        Temperature = [float(input("Enter the temperature (T): "))]  # If single temperature, prompt user for that temperature
    else:
        start_temp = float(input("Enter the starting temperature: "))  # If multiple temperatures, prompt user for temperature range
        stop_temp = float(input("Enter the stopping temperature: "))
        step_temp = float(input("Enter the temperature step size: "))
        Temperature = np.arange(start_temp, stop_temp + step_temp, step_temp)  # Generate array of temperatures
    return B, N, n, Temperature

# Monte Carlo sweep implementation
def metropolis_algorithm(lattice, Temp, n, N, B):   
    for m in range(n):
        i = random.randrange(N)  # Choose random row
        j = random.randrange(N)  # Choose random column
        energy_diff = delta_energy(lattice, i, j, N, B)  # Calculate change in energy
        if energy_diff <= 0:  # If change in energy is negative, accept move and flip spin
            lattice[i, j] = -lattice[i, j]
        elif random.random() < np.exp(-energy_diff/(boltz_cont * Temp)):  # If change in energy is positive, accept move with certain probability
            lattice[i, j] = -lattice[i, j]
    return lattice

# Compute physical quantities
def calc_parameters(lattice, T, n, N, B):
    magnetism=0
    energy=0
    magnetism_squared=0
    energy_squared=0
    for p in range(n):
        lattice = metropolis_algorithm(lattice, T, 1, N, B)
        E = calc_energy(lattice, N, B)
        M = calc_mag(lattice, N)
        energy += E
        magnetism += M
        energy_squared += E*E
        magnetism_squared += M*M
    average_energy = energy / n
    mag = magnetism / n
    heat_cap = (energy_squared / n - (energy / n) ** 2) / (boltz_cont*T**2)  # Specific heat capacity
    susceptibility = (magnetism_squared / n - (magnetism / n) ** 2) / (boltz_cont*T)  # Susceptibility
    return average_energy, mag, heat_cap, susceptibility

def run_simulation(B, N, n, Temperature):
    lattice = np.random.choice([1,-1], size=(N, N))  # Initialize lattice with random spins
    
    # Initialize arrays to store results
    mag = np.zeros(len(Temperature))
    average_energy = np.zeros(len(Temperature))
    CV = np.zeros(len(Temperature))
    susceptibility_arr = np.zeros(len(Temperature))

    # Open a text file to store the data
    with open("simulation_results.txt", "w") as file:
        file.write("Temperature\tAverage Energy\tMagnetization\tSpecific Heat\tSusceptibility\n")
        
        # Print simulation details to the file
        file.write(f"Lattice Size: {N}x{N}\n")
        file.write(f"External Magnetic Field (B): {B}\n")
        file.write(f"Metropolis Step: {n}\n")

        start = time.time()  # Start timing

        # Simulate at particular temperatures (T) and compute physical quantities (Energy, heat capacity, magnetization, susceptibility)
        for ind, T in enumerate(tqdm(Temperature)):
            lattice = metropolis_algorithm(lattice, T, n, N, B)  # Perform Monte Carlo sweeps
            average_energy[ind], mag[ind], CV[ind], susceptibility_arr[ind] = calc_parameters(lattice, T, n, N, B)  # Compute physical quantities
            
            # Write the results to the file
            file.write(f"{T}\t{average_energy[ind]}\t{mag[ind]}\t{CV[ind]}\t{susceptibility_arr[ind]}\n")

        end = time.time()
        time_taken = (end - start)
        file.write(f"Time taken: {time_taken} seconds\n")

    # Plot results
    plot_results(Temperature, average_energy, mag, CV, susceptibility_arr, B, N, n)




# Function to display menu options
def display_menu():
    print("\n===== 2D Ising Model Simulation Menu =====")
    print("1. Run Simulation")
    print("2. Exit")

# Main function to run the program
def main():
    while True:
        display_menu()
        choice = input("Enter your choice: ")
        
        if choice == '1':
            B, N, n, Temperature = get_parameters()
            run_simulation(B, N, n, Temperature)
        elif choice == '2':
            print("Exiting the program. Goodbye!")
            break
        else:
            print("Invalid choice. Please enter a valid option.")

if __name__ == "__main__":
    main()
