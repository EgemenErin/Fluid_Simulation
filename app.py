import pygame
import numpy as np

N = 64  #Number of cells in each dimension (simulation grid cells are 1...N)
SIZE = N + 2  # Add boundaries on each side
DT = 0.1  #Time step
DIFF = 0.0  #Diffusion rate (set to 0 for no density diffusion)
VISC = 0.0  #Viscosity (set to 0 for inviscid flow)
ITER = 4  # Number of iterations for the linear solver
dens = np.zeros((SIZE, SIZE), dtype=float)  #Density
dens_prev = np.zeros((SIZE, SIZE), dtype=float)  #Density source
u = np.zeros((SIZE, SIZE), dtype=float)  #X-velocity
v = np.zeros((SIZE, SIZE), dtype=float)  #Y-velocity
u_prev = np.zeros((SIZE, SIZE), dtype=float)  #Velocity source in x
v_prev = np.zeros((SIZE, SIZE), dtype=float)  # elocity source in y

def set_bnd(b, x):
    """ Sets the boundary conditions. For b==1, we are dealing with horizontal velocities;for b==2, vertical; and for b==0, the scalar (density). """
    for i in range(1, N + 1):
        if b == 1:
            x[0, i] = -x[1, i]
            x[N + 1, i] = -x[N, i]
        else:
            x[0, i] = x[1, i]
            x[N + 1, i] = x[N, i]
        if b == 2:
            x[i, 0] = -x[i, 1]
            x[i, N + 1] = -x[i, N]
        else:
            x[i, 0] = x[i, 1]
            x[i, N + 1] = x[i, N]
    x[0, 0] = 0.5 * (x[1, 0] + x[0, 1])
    x[0, N + 1] = 0.5 * (x[1, N + 1] + x[0, N])
    x[N + 1, 0] = 0.5 * (x[N, 0] + x[N + 1, 1])
    x[N + 1, N + 1] = 0.5 * (x[N, N + 1] + x[N + 1, N])


def add_source(x, s, dt):
    """adds the source term s into x over a timestep dt."""
    x += dt * s


def diffuse(b, x, x0, diff, dt):
    """Diffuses quantity x0 into x."""
    a = dt * diff * N * N
    for k in range(ITER):
        for i in range(1, N + 1):
            for j in range(1, N + 1):
                x[i, j] = (x0[i, j] + a * (x[i - 1, j] + x[i + 1, j] +
                                           x[i, j - 1] + x[i, j + 1])) / (1 + 4 * a)
        set_bnd(b, x)


def advect(b, d, d0, u, v, dt):
    """Transports (advects) quantity d0 into d using thee velocity field (u,v). It uses a backward-tracing method."""
    dt0 = dt * N
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            x = i - dt0 * u[i, j]
            y = j - dt0 * v[i, j]
            if x < 0.5: x = 0.5
            if x > N + 0.5: x = N + 0.5
            if y < 0.5: y = 0.5
            if y > N + 0.5: y = N + 0.5
            i0 = int(x)
            i1 = i0 + 1
            j0 = int(y)
            j1 = j0 + 1
            s1 = x - i0
            s0 = 1 - s1
            t1 = y - j0
            t0 = 1 - t1
            d[i, j] = (s0 * (t0 * d0[i0, j0] + t1 * d0[i0, j1]) +
                       s1 * (t0 * d0[i1, j0] + t1 * d0[i1, j1]))
    set_bnd(b, d)


def project(u, v, p, div):
    h = 1.0 / N
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            div[i, j] = -0.5 * h * (u[i + 1, j] - u[i - 1, j] +
                                    v[i, j + 1] - v[i, j - 1])
            p[i, j] = 0
    set_bnd(0, div)
    set_bnd(0, p)
    for k in range(ITER):
        for i in range(1, N + 1):
            for j in range(1, N + 1):
                p[i, j] = (div[i, j] + p[i - 1, j] + p[i + 1, j] +
                           p[i, j - 1] + p[i, j + 1]) / 4
        set_bnd(0, p)
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            u[i, j] -= 0.5 * (p[i + 1, j] - p[i - 1, j]) / h
            v[i, j] -= 0.5 * (p[i, j + 1] - p[i, j - 1]) / h
    set_bnd(1, u)
    set_bnd(2, v)


def vel_step(u, v, u0, v0, visc, dt):
    add_source(u, u0, dt)
    add_source(v, v0, dt)
    u0[:, :] = u.copy()
    v0[:, :] = v.copy()
    diffuse(1, u, u0, visc, dt)
    diffuse(2, v, v0, visc, dt)
    p = np.zeros_like(u)
    div = np.zeros_like(u)
    project(u, v, p, div)
    u0[:, :] = u.copy()
    v0[:, :] = v.copy()
    advect(1, u, u0, u0, v0, dt)
    advect(2, v, v0, u0, v0, dt)
    project(u, v, p, div)


def dens_step(x, x0, u, v, diff, dt):
    add_source(x, x0, dt)
    x0[:, :] = x.copy()
    diffuse(0, x, x0, diff, dt)
    x0[:, :] = x.copy()
    advect(0, x, x0, u, v, dt)




def main():
    pygame.init()
    window_size = 512
    screen = pygame.display.set_mode((window_size, window_size))
    pygame.display.set_caption("Interactive Fluid Simulation")
    clock = pygame.time.Clock()

    cell_size = window_size // N

    running = True
    last_mouse = None

    global dens, dens_prev, u, u_prev, v, v_prev

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
        mouse_pressed = pygame.mouse.get_pressed()
        mouse_pos = pygame.mouse.get_pos()
        grid_x = int(mouse_pos[0] / cell_size) + 1
        grid_y = int(mouse_pos[1] / cell_size) + 1

        if mouse_pressed[0]:
            if last_mouse is not None:
                dx = mouse_pos[0] - last_mouse[0]
                dy = mouse_pos[1] - last_mouse[1]
            else:
                dx, dy = 0, 0

            #can adjust these numbers to change the strength of the input to play around.
            force = 100.0
            u_prev[grid_x, grid_y] = dx * force
            v_prev[grid_x, grid_y] = dy * force
            dens_prev[grid_x, grid_y] = 100.0
        last_mouse = mouse_pos
        vel_step(u, v, u_prev, v_prev, VISC, DT)
        dens_step(dens, dens_prev, u, v, DIFF, DT)
        u_prev.fill(0)
        v_prev.fill(0)
        dens_prev.fill(0)

        screen.fill((0, 0, 0))
        for i in range(1, N + 1):
            for j in range(1, N + 1):
                d = dens[i, j]
                if d > 255:
                    d = 255
                color = (int(d), int(d), int(d))
                rect = pygame.Rect((i - 1) * cell_size, (j - 1) * cell_size, cell_size, cell_size)
                pygame.draw.rect(screen, color, rect)

        pygame.display.flip()
        clock.tick(60)

    pygame.quit()


if __name__ == '__main__':
    main()
