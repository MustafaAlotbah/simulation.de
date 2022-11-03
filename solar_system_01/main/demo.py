"""
    Demo on the simulation of solar systems for simulation.de visualizations
    Mustafa Alotbah 2022
"""


from solar_system_01.simulation import Simulation, Particle
from solar_system_01.visualization import Visualization

if __name__ == "__main__":
    print("Simulation starting...")

    sim = Simulation()

    sim.particles.append(Particle(
        label="Sun", color="#dcd427", mass=100_000.0, x=0.0, y=0.0
    ))

    sim.particles.append(Particle(
        label="Earth", color="#0092cc", mass=1_000.0, x=0.0, y=2.0, vx=0.6
    ))

    sim.particles.append(Particle(
        label="Moon", color="white", mass=351.0, x=0.0, y=2.4, vx=0.632
    ))

    sim.particles.append(Particle(
        label="Asteroid", color="#ff3333", mass=0.005, x=0.0, y=6.55, vx=0.152
    ))

    ani = Visualization(
        sim,
        bounds=[-5, 5, -12, 8],
        simulation_steps=3001,
        background_color='#020203',
        foreground_color='white',
    )

    ani.show()

    # import os
    # num_files = len([name for name in os.listdir('../output') if os.path.isfile('../output/' + name)])
    # ani.save_video(f'../output/{num_files}.mp4')

    print("Done.")
