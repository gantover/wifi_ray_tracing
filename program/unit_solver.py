from utils import *
from wall import Wall  # importing static methods for coefficients calculation
from data import P_MAX_CL


@ti.func
def intersect(u, n, p1, q1, p2, q2):
    # checks if there is an intersection between wall (u,n,p1,q1)
    # and a ray (p2, q2)
    value = False
    if tm.sign((q2 - p1).dot(n)) != tm.sign((p2 - p1).dot(n)):
        t = find_intersection(p1, u, p2, q2)
        ip = Wall.point_on_wall(p1, u, t)
        if tm.sign((ip - p1).dot(u)) != tm.sign((ip - q1).dot(u)):
            value = True
    return value


@ti.func
def dist(r0, u, p):
    return ti.abs(u[1] * p[0] - u[0] * p[1] - u[1] * r0[0] + u[0] * r0[1])


@ti.func
def get_next_tx(r0, u, n, tx):
    # gets the symmetric of tx by the wall plane (r0, u, n)
    return tx - 2 * n * tm.sign((tx - r0).dot(n)) * dist(r0, u, tx)


@ti.func
def find_intersection(r0, u, p2, q2):
    # this finds t, which will give the intersection point of the ray p2,q2
    # and the wall (r0,u) by intersection = t * u
    l = q2 - p2
    dx = l[0]
    dy = l[1]
    t = (dy * (r0[0] - q2[0]) - dx * (r0[1] - q2[1])) / (dx * u[1] - dy * u[0])
    return t


@ti.func
def bounce_cond(r0, n, p2, q2):
    # basic condition to know if the bounce is physically possible
    return tm.sign(n.dot(p2 - r0)) == tm.sign(n.dot(q2 - r0))


@ti.func
def wall_transmission(world, p2, q2, index1=-1, index2=-1):
    # see which walls are in the way of the ray defined by p2, q2
    # and calculates the total transmission coefficient for this ray
    # index1, index2 are walls that we don't want to take into account
    if index1 > index2:
        index_temp = index1
        index1 = index2
        index2 = index_temp

    transmission_factor_msq = ti.f32(1.0)
    ray_normal = (q2 - p2).normalized()
    for i in range(0, index1):
        p1, q1 = world.r0[i], world.r1[i]
        u, n = world.u[i], world.n[i]
        if intersect(u, n, p1, q1, p2, q2):
            transmission_factor_msq *= Wall.get_tn2(world, i, ray_normal)
    for i in range(index1 + 1, index2):
        p1, q1 = world.r0[i], world.r1[i]
        u, n = world.u[i], world.n[i]
        if intersect(u, n, p1, q1, p2, q2):
            transmission_factor_msq *= Wall.get_tn2(world, i, ray_normal)
    for i in range(index2 + 1, world.m):
        p1, q1 = world.r0[i], world.r1[i]
        u, n = world.u[i], world.n[i]
        if intersect(u, n, p1, q1, p2, q2):
            transmission_factor_msq *= Wall.get_tn2(world, i, ray_normal)
    return transmission_factor_msq


@ti.func
def calculate_power(world, tx, rx):
    # we avoid the problem of not being in "distant fields" hypothesis
    # by limiting the power value to its base value at a distance of 1m
    prx_temp = tm.clamp(1.0 / (tm.pow((rx - tx).norm(), 2)), 0.0, P_MAX_CL)

    trans0_0 = wall_transmission(world, rx, tx)
    prx = prx_temp * trans0_0

    for i in range(world.m):
        transmission_factor_msq = ti.f32(1.0)
        reflexion_factor_msq = ti.f32(1.0)
        r0_i, r1_i = world.r0[i], world.r1[i]
        u_i, n_i = world.u[i], world.n[i]
        tx1 = get_next_tx(r0_i, u_i, n_i, tx)
        if bounce_cond(r0_i, n_i, tx, rx):
            t = find_intersection(r0_i, u_i, tx1, rx)
            ip = Wall.point_on_wall(r0_i, u_i, t)
            if tm.sign((ip - r0_i).dot(u_i)) != tm.sign((ip - r1_i).dot(u_i)):
                reflexion_factor_msq *= Wall.get_rn2(world, i, (ip - tx).normalized())
                transmission_factor_msq *= wall_transmission(world, tx, ip, i) \
                                           * wall_transmission(world, ip, rx, i)
                prx_temp = tm.clamp(1.0 / ((rx - tx1).norm() ** 2), 0.0, P_MAX_CL)
                prx += prx_temp * transmission_factor_msq * reflexion_factor_msq

        for j in range(world.m):
            if i == j:
                continue
            transmission_factor_msq = ti.f32(1.0)
            reflexion_factor_msq = ti.f32(1.0)
            r0_j, r1_j = world.r0[j], world.r1[j]
            u_j, n_j = world.u[j], world.n[j]
            if not bounce_cond(r0_j, n_j, tx1, rx):
                continue
            tx2 = get_next_tx(r0_j, u_j, n_j, tx1)
            t2 = find_intersection(r0_j, u_j, tx2, rx)
            ip2 = Wall.point_on_wall(r0_j, u_j, t2)
            if tm.sign((ip2 - r0_j).dot(u_j)) == tm.sign((ip2 - r1_j).dot(u_j)):
                continue
            if not bounce_cond(r0_i, n_i, tx, ip2):
                continue
            t1 = find_intersection(r0_i, u_i, tx1, ip2)
            ip1 = Wall.point_on_wall(r0_i, u_i, t1)
            if tm.sign((ip1 - r0_i).dot(u_i)) == tm.sign((ip1 - r1_i).dot(u_i)):
                continue
            reflexion_factor_msq *= Wall.get_rn2(world, i, (ip1 - tx).normalized()) \
                                    * Wall.get_rn2(world, j, (ip2 - ip1).normalized())
            transmission_factor_msq *= wall_transmission(world, tx, ip1, i) \
                                       * wall_transmission(world, ip1, ip2, j, i) \
                                       * wall_transmission(world, ip2, rx, j)
            prx_temp = tm.clamp(1.0 / ((rx - tx2).norm() ** 2), 0.0, P_MAX_CL)
            prx += prx_temp * reflexion_factor_msq * transmission_factor_msq
    return prx
