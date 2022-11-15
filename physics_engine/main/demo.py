"""
    Demo on the simulation of spring-mass system for simulation.de visualizations
    Mustafa Alotbah 2022
"""


from physics_engine.simulation import Simulation, Particle, Spring
from physics_engine.visualization import Visualization

if __name__ == "__main__":
    print("Simulation starting...")

    sim = Simulation(dt=0.005)

    # first pendulum --------------------
    # Anchor
    sim.particles.append(Particle(
        label="Anchor", color="lightgrey", mass=0.5, x=0.0, y=0.0, is_fixed=True
    ))
    sim.particles[-1].apply_earth_gravity = False

    # -- First Ball
    sim.particles.append(Particle(
        label="Ball 1", color="#0092cc", mass=0.5, x=-1.5, y=-0.5, vx=0.0, is_fixed=False
    ))
    sim.springs.append(Spring(sim.particles[0], sim.particles[1]))
    sim.springs[-1].rest_length = 1.5811388300841898

    # -- Second Ball
    sim.particles.append(Particle(
        label="Ball 2", color="#dcd427", mass=0.5, x=-1.5, y=-0.5+1, vx=0.0, is_fixed=False
    ))
    sim.springs.append(Spring(sim.particles[1], sim.particles[2]))

    # # Second pendulum --------------------
    # # same anchor
    #
    # # -- First Ball
    # sim.particles.append(Particle(
    #     label="Ball 1\'", color="#ff3333", mass=1.5, x=-1.5, y=-0.5, vx=0.0, is_fixed=False
    # ))
    # sim.springs.append(Spring(sim.particles[0], sim.particles[3]))
    # sim.springs[-1].rest_length = 1.5811388300841898
    #
    # # -- Second Ball
    # sim.particles.append(Particle(
    #     label="Ball 2\'", color="#88bb44", mass=0.5, x=-1.5, y=-0.5+1, vx=0.0, is_fixed=False
    # ))
    # sim.springs.append(Spring(sim.particles[3], sim.particles[4]))

    # end --------------------------------
    ani = Visualization(
        sim,
        bounds=[-3, 3, -10, 2],
        simulation_steps=3001,
        background_color='#020203',
        foreground_color='white',
    )

    ani.show()

    # import os
    # num_files = len([name for name in os.listdir('../output') if os.path.isfile('../output/' + name)])
    # ani.save_video(f'../output/{num_files}.mp4')

    print("Done.")
