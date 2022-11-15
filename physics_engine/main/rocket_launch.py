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

    sim = Simulation(dt=0.003)

    # Rocket Simulation Parameters
    engine_mass = 0.02
    rocket_mass = 0.05
    thrust_force_before_tilt = 9.81 * 0.07          # these forces happen per particle emission
    thrust_force_after_tilt = 9.81 * 0.17
    time_to_tilt = 0.6
    tilt_duration = 0.1
    tilt_in_degrees = 5
    particle_emission_per_step = 2
    particle_age_decrease_per_step = 0.005

    # <editor-fold desc="The Design of The Rocket">
    # The first point is also the engine
    sim.particles.append(Particle(
        label="Engine", color="#ff3333", mass=engine_mass, x=0.0, y=0.0, has_legend=True, trace_alpha=0.6, radius=0.05,
        display_forces=False
    ))
    sim.particles.append(Particle(
        label="Rocket", color="lightgrey", mass=rocket_mass/5, x=-0.25, y=0.0, has_legend=True, trace_alpha=0.02,
        display_forces=False
    ))
    sim.particles.append(Particle(
        label="Anchor 1", color="lightgrey", mass=rocket_mass/5, x=0.25, y=0.0, has_legend=False, trace_alpha=0.02,
        display_forces=False
    ))
    sim.particles.append(Particle(
        label="Anchor 1", color="lightgrey", mass=rocket_mass/5, x=-0.25, y=1.0, has_legend=False, trace_alpha=0.02,
        display_forces=False
    ))
    sim.particles.append(Particle(
        label="Anchor 1", color="lightgrey", mass=rocket_mass/5, x=0.25, y=1.0, has_legend=False, trace_alpha=0.02,
        display_forces=False
    ))
    sim.particles.append(Particle(
        label="Anchor 1", color="lightgrey", mass=rocket_mass/5, x=0.0, y=1.25, has_legend=False, trace_alpha=0.02,
        display_forces=False
    ))

    # Rocket edges
    sim.springs.append(Spring(sim.particles[0], sim.particles[1], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[0], sim.particles[2], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[3], sim.particles[1], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[4], sim.particles[2], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[3], sim.particles[2], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[4], sim.particles[1], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[4], sim.particles[3], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[4], sim.particles[5], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[3], sim.particles[5], stiffness=50.0, damping=1.0))
    sim.springs.append(Spring(sim.particles[0], sim.particles[5], stiffness=50.0, damping=1.0))
    # </editor-fold>

    # <editor-fold desc="The ground">
    sim.particles.append(Particle(
        label="Ground", color="#88bb44", mass=0.1, x=1.0, y=-0.25, has_legend=True, trace_alpha=0.02, radius=0.07, is_fixed=True,
        display_forces=False
    ))
    sim.particles[-1].apply_earth_gravity = False
    sim.particles.append(Particle(
        label="", color="#88bb44", mass=0.1, x=-1.0, y=-0.25, has_legend=False, trace_alpha=0.02, radius=0.07, is_fixed=True,
        display_forces=False
    ))
    sim.particles[-1].apply_earth_gravity = False
    sim.springs.append(Spring(sim.particles[-1], sim.particles[-2], stiffness=50.0))
    # </editor-fold>

    # <editor-fold desc="The Particles emitted from the engine">
    # a color map to change the particles' colors by age alter
    norm = colors.TwoSlopeNorm(vcenter=0.5, vmin=0.0, vmax=1.0)
    sim.color_map = cm.ScalarMappable(cmap='inferno', norm=norm)

    # add dead particles (to allocate their memory [see Matplotlib Animation])
    particles_index = len(sim.particles)
    for i in range(200):
        sim.particles.append(Particle(
            label="Anchor 1", color="lightgrey", mass=0.02, x=0.5, y=0.0, has_legend=False, trace_alpha=0.02, display_forces=False
        ))
        sim.particles[-1].age = 0
        sim.particles[-1].alpha = 0
        sim.particles[-1].apply_earth_gravity = False  # This is only temporary

    # </editor-fold>

    # <editor-fold desc="Rocket Simulation Functionality">
    # This is the main particles emitter functionality
    sim.rocket_orientation = 0

    # bring dead particles to live one after another, and make the others older
    def callback(sim):
        # dead particles to life every other step
        particles_step = 0
        if (sim.time//sim.dt) % particle_emission_per_step < 0.001:
            for particle in sim.particles[particles_index:]:
                if particle.age <= 0:
                    # thrust magnitude and orientation
                    thrust_magnitude = thrust_force_before_tilt
                    if time_to_tilt + tilt_duration <= sim.time:
                        thrust_magnitude = thrust_force_after_tilt

                    theta = sim.rocket_orientation * np.pi / 180
                    if time_to_tilt <= sim.time <= time_to_tilt + tilt_duration:
                        theta -= tilt_in_degrees * np.pi / 180

                    thrust_y = thrust_magnitude * np.cos(theta)
                    thrust_x = thrust_magnitude * np.sin(theta)

                    # emit from the engine's point
                    particle.position = sim.particles[0].position.copy()

                    # direction with slight randomness (due to air pressure but not simulated here)
                    particle.velocity = np.array([(np.random.rand(1)[0]-0.5)*4.5 - thrust_x, -np.random.rand(1)[0]*2.0 - thrust_y])
                    particle.positions = [tuple(particle.position)]
                    particle.alpha = 1.0
                    particle.age = 1.0
                    sim.particles[-1].apply_earth_gravity = True
                    particles_step += 1

                    # add thrust to the rocket's engine
                    sim.particles[0].forces.append(np.array([thrust_x, thrust_y]))

                if particles_step >= 2:  # one particle at a time
                    break

        # update the ages and the colors accordingly
        for particle in sim.particles[particles_index:]:
            if particle.age > 0.00:
                particle.age -= particle_age_decrease_per_step
                particle.alpha = min(max(0.0, particle.age), 1.0)
                particle.color = sim.color_map.to_rgba(0.1 + 0.75 * particle.age)
                particle.radius = 0.03 + 0.05 * (1-particle.age)        # radius gets larger by decreasing age

        # collision with ground
        for particle in sim.particles:
            if particle.position[1] < -0.25 and -1 < particle.position[0] < 1:
                particle.velocity[1] = -particle.velocity[1]    # reflects its velocity on the y-axis
                if 'age' in dir(particle):
                    particle.age -= 0.1                # a collision makes the particle age more instantly
                    particle.velocity[0] *= 2          # also diffuses its orientation

    sim.callback = callback
    # </editor-fold>


    def stats_callback(sim):
        a = sim.particles[0].position
        b = sim.particles[5].position
        x = (a-b)/np.linalg.norm(a-b)
        sim.rocket_orientation = np.arctan(x[0]/x[1]) * 180 / np.pi
        return (
            f"Engine's velocity: {np.linalg.norm(sim.particles[0].velocity):.2f} m/s\n"
            f"Tilt: {sim.rocket_orientation:.2f}°"
        )

    sim.stats_callback = stats_callback


    # end --------------------------------
    ani = Visualization(
        sim,
        bounds=[-4, 4, -3, 13],
        simulation_steps=1501,
        background_color='#020203',
        foreground_color='white',
        frames_skipped=0
    )

    ani.show()


    # import os
    # path = "physics_engine/output"
    # num_files = len([name for name in os.listdir(path) if os.path.isfile(path + '/' + name)])
    # print(num_files)
    # ani.save_video(f'{path}/{num_files}.mp4')
    #
    print("Done.")
