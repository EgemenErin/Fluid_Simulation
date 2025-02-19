# Interactive Fluid Simulation

This project is an interactive fluid simulation implemented in Python. The simulation uses a grid-based approach inspired by Jos Stam’s "Stable Fluids" algorithm to model the motion of fluids in a two-dimensional space.
I cannot get credit for this project hence i followed a youtube tutorial to make it.

## Features

- **Fluid Dynamics:** Simulates fluid behavior such as diffusion, advection, and incompressibility using a simplified version of the Navier–Stokes equations.
- **User Interaction:** Click and drag the mouse to inject density and velocity, creating dynamic fluid motion.
- **Visualization:** The simulation renders a colorful fluid (blue on a white background) that reacts in real time to user input.

## How It Works

1. **Grid Setup:**  
   The fluid is represented on a 2D grid. Each cell in the grid stores properties like density and velocity.

2. **Physics Processes:**  
   - **Add Source:** External forces (like mouse input) add velocity and density.
   - **Diffuse:** Simulates viscosity and the natural spreading of density.
   - **Advect:** Moves the density and velocity along the flow of the fluid.
   - **Project:** Enforces incompressibility, ensuring the fluid remains continuous with no gaps or piles.

3. **Rendering:**  
   The updated density field is rendered in real time, displaying the fluid's motion as it evolves.

## Getting Started

### Prerequisites

- Python 3.x
- [NumPy](https://numpy.org/)
- [Pygame](https://www.pygame.org/)

### Installation

**Clone the repository:**
bash git clone https://github.com/EgemenErin/Fluid_Simulation
 
**Navigate to the project directory:**
cd Fluid_Simulation

**Install the required packages:**
pip install numpy pygame

**Run the simulation by executing:**
python fluid_simulation.py

### Acknowledgments

I would like to extend my sincere thanks to my mathematician friend for their invaluable help in understanding the mathematical foundations behind fluid dynamics and the Navier–Stokes equations. 
Their insights played a significant role in understanding the project. I wouldnt be able to follow and figure out the equations if it werent him. He helped me lots of amounts over these past three months with these projects and even contributed while i was dealing with my finals.

