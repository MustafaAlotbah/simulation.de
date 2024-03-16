"""
    Demo on the simulation of spring-mass system for simulation.de visualizations
    Mustafa Alotbah 2022
"""


from physics_engine.simulation import Simulation, Particle, Spring, INTEGRATOR_LEAPFROG, INTEGRATOR_RK4
from physics_engine.visualization import Visualization

if __name__ == "__main__":
    print("Simulation starting...")

    sim = Simulation(dt=0.005, integrator=INTEGRATOR_RK4)
    x_increment = 0.075
    length = 3.5
    num_joints = int(length / x_increment)

    # first pendulum --------------------
    # Anchor
    placement_x = 0.0
    sim.particles.append(Particle(
        label="Anchor", color="#0092cc", mass=0.5, radius=0.1, x=placement_x, y=0.0, is_anchored=True, display_forces=False, trace_alpha=0
    ))
    sim.particles[-1].apply_earth_gravity = False

    # -- add connections
    for i in range(num_joints):
        placement_x += x_increment
        sim.particles.append(Particle(
            has_legend=False, color="lightgrey", mass=0.05, x=placement_x, y=0.0, vx=0.0, is_anchored=False, display_forces=False, trace_alpha=0
        ))
        sim.springs.append(Spring(sim.particles[-2], sim.particles[-1], dt=sim.dt))

    # -- Final Ball
    placement_x += x_increment
    sim.particles.append(Particle(
        label="Mass", color="#dcd427", mass=1.0, radius=0.1, x=placement_x, y=-0.0, vx=0.0, is_anchored=False, display_forces=False
    ))
    sim.springs.append(Spring(sim.particles[-2], sim.particles[-1], dt=sim.dt))

    # end --------------------------------
    ani = Visualization(
        sim,
        bounds=[-4, 4, -11, 3],
        simulation_steps=3001,
        background_color='#020203',
        foreground_color='white',
        frames_skipped=1
    )

    ani.show()

    # import os
    # num_files = len([name for name in os.listdir('../output') if os.path.isfile('../output/' + name)])
    # ani.save_video(f'../output/{num_files}.mp4')

    print("Done.")
