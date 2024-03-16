"""
    Demo on the simulation of spring-mass system for simulation.de visualizations
    Mustafa Alotbah 2022
"""

import numpy as np
from physics_engine.simulation import Simulation, Particle, Spring
from physics_engine.visualization import Visualization
from matplotlib import colors, cm


if __name__ == "__main__":
    print("Simulation starting...")

    sim = Simulation(dt=0.005)

    mass = 0.01

    # a color map to change the particles' colors by age alter
    norm = colors.TwoSlopeNorm(vcenter=0.5, vmin=0.0, vmax=1.0)
    sim.color_map = cm.ScalarMappable(cmap='inferno', norm=norm)

    #
    # add dead particles (to allocate their memory [see Matplotlib Animation])
    for i in range(30):
        sim.particles.append(Particle(
            label="Anchor 1", color="lightgrey", mass=mass, x=0.5, y=0.0, has_legend=False, trace_alpha=0.02, display_forces=False
        ))
        sim.particles[-1].age = 0
        sim.particles[-1].alpha = 0
        sim.particles[-1].apply_earth_gravity = False  # This is only temporary

    # This is the main particles emitter functionality
    # bring dead particles to live one after another, and make the others older
    def callback(sim_state):
        # dead particles to life every half second
        particles_step = 0
        if (sim_state.uptime//sim_state.dt) % 2 < 0.001:
            for particle in sim_state.particles:
                if particle.age <= 0:
                    # The emitter is the third particle
                    particle.position = np.array([0.0, 0.0])
                    particle.velocity = np.array([(np.random.rand(1)[0]-0.5)*3.0, np.random.rand(1)[0]*2.0 + 4.0])

                    particle.positions = [tuple(particle.position)]
                    particle.alpha = 1.0
                    particle.age = 1.0
                    sim_state.particles[-1].apply_earth_gravity = True
                    particles_step += 1

                if particles_step >= 2:
                    break

        # update the ages and the colors accordingly
        for particle in sim_state.particles:
            if particle.age > 0.00:
                particle.age -= 0.005
                particle.alpha = min(max(0.0, particle.age), 1.0)
                particle.color = sim.color_map.to_rgba(0.1 + 0.75 * particle.age)
                particle.radius = 0.03 + 0.04 * (1-particle.age)



    sim.callback = callback


    # end --------------------------------
    ani = Visualization(
        sim,
        bounds=[-3, 3, -6, 6],
        simulation_steps=3501,
        background_color='#020203',
        foreground_color='white',
        frames_skipped=0
    )

    ani.show()

    # import os
    # num_files = len([name for name in os.listdir('../output') if os.path.isfile('../output/' + name)])
    # ani.save_video(f'../output/{num_files}.mp4')

    print("Done.")
