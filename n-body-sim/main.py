import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation

# Important variables declared
# l - determines number of time steps
# num_time_steps - number of discrete time segments
# t_max - the length of simulation
# dt - delta of time per time segment
l = 10
num_time_steps = 2**l + 1
t_max = 50
dt = t_max/(num_time_steps-1)

# Set number of particles in simulation
num_particles = 8

# Indices of array
x = 0
y = 1
z = 2

# Initialize Masses
m = np.random.randint(15, 50, size=num_particles)

# Attach 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

plot_lim = 10


def set_fig():
    ax.cla()
    # Set the axes properties
    ax.set_xlim3d([-1*plot_lim, plot_lim])
    ax.set_xlabel('X')

    ax.set_ylim3d([-1*plot_lim, plot_lim])
    ax.set_ylabel('Y')

    ax.set_zlim3d([-1*plot_lim, plot_lim])
    ax.set_zlabel('Z')

    ax.set_title('3D Test')


# Hacky way of updating lines by clearing graph and re plotting
def update_lines(t, r, lines):
    set_fig()
    ax.set_title('3D Test')
    for i in range(0, num_particles):
        ax.plot(r[i, x, t], r[i, y, t], r[i, z, t], 'o')

    return lines


def randomize_initial_conditions(r, v):
    for i in range(0, num_particles):
        r[i, :, 0] = (np.random.rand(3)*2-1)*8
        r[i, :, 1] = r[i, :, 0]
        v[i] = (np.random.rand(3)*2 - 1)*3


def main():
    # Initialize displacement and velocity arrays
    r = np.zeros((num_particles, 3, num_time_steps))
    v = np.zeros((num_particles, 3))

    randomize_initial_conditions(r, v)

    # Calculate position for t = 2
    for i in range(0, num_particles):
        calc_position_for_t2(r, i, m, v)

    # Calculate position for all t
    for t in range(3, num_time_steps):
        for i in range(0, num_particles):
            for k in range(0, 3):
                calculate_for_t(r, k, i, m, t)

    # Initialize lines for each particle
    lines = [ax.plot(r[0, x, 0], r[0, y, 0], r[0, z, 0], 'o')]

    num_iterations = 500
    # Creating the Animation object
    ani = animation.FuncAnimation(
        fig, update_lines, num_iterations, fargs=(r, lines), interval=50)

    # Uncomment to save as gif
    # ani.save("test.gif")
    plt.show()


def calc_position_for_t2(r, i, m, v):
    t = 2
    next_x = r[i, x, 0] + dt*v[i, x] + \
        (1/2)*(dt**2)*sum_for_dirn(r, x, i, m, 0)
    next_y = r[i, y, 0] + dt*v[i, y] + \
        (1/2)*(dt**2)*sum_for_dirn(r, y, i, m, 0)
    next_z = r[i, z, 0] + dt*v[i, z] + \
        (1/2)*(dt**2)*sum_for_dirn(r, z, i, m, 0)

    r[i, x, t] = next_x
    r[i, y, t] = next_y
    r[i, z, t] = next_z


# r - displacement vector. 3d matrix where each row is a particle, each column is for x y z and each layer is for each time segment
# dirn - direction calculated. e.g. x y or z
# i - particle i
# m - masses of particles
# t - time slice
def calculate_for_t(r, dirn, i, m, t):
    next_position = 2*r[i, dirn, t-1] - r[i, dirn, t-2] + \
        (dt**2)*sum_for_dirn(r, dirn, i, m, t-1)
    r[i, dirn, t] = next_position


def sum_for_dirn(r, dirn, i, m, t):
    sum = 0
    for j in range(0, num_particles):
        if (i != j):
            sum = sum + ((m[j]*(r[j, dirn, t] - r[i, dirn, t])) /
                         (get_rij_for_time(r, i, j, t)))

    return sum


def get_rij_for_time(r, i, j, t):
    xdiff = np.power(r[j, x, t] - r[i, x, t], 2)
    ydiff = np.power(r[j, y, t] - r[i, y, t], 2)
    zdiff = np.power(r[j, z, t] - r[i, z, t], 2)

    return np.power(xdiff + ydiff + zdiff, 3/2)


if __name__ == "__main__":
    main()
