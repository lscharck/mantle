#include "mantle.h"
extern "C" {
	#include "rebound.h"
}

template<typename T, typename F>
void liftOptional(std::optional<T> val, F f) {
    if (val.has_value()) {f(val.value());}
}

Mantle::Mantle(m_sim_params user_params, userParticles userFn) {

    reb_sim_params = user_params;
    userParticlesFn = userFn;
    has_additional_forces = false;
    has_ode = false;

}

Mantle::Mantle(m_sim_params user_params, std::vector<Particle> userParticles) {

    reb_sim_params = user_params;
    has_additional_forces = false;
    has_ode = false;
    particles = userParticles;

}

Mantle::~Mantle() {

    reb_simulation_free(r);

}

void Mantle::setAdditionalForce(void) {

    has_additional_forces = true;

}

void Mantle::setCollisionContext(m_coll_context userCollisionContext) {

    collisionContext = userCollisionContext;

}

void Mantle::setOde(void) {

    has_ode = true;

}

long Mantle::getCollisionLog() const {

    return r->collisions_log_n;

}

void Mantle::initMantle(void) {

    r = reb_simulation_create();
    updateParams();
    if (userParticlesFn.has_value()) {auto f = userParticlesFn.value(); particles = f();}
    addHash();
    if (collisionContext.has_value()) {addCollision(collisionContext.value());};
    if (has_additional_forces) {addAdditionalForce();}
    if (has_ode) {addOde();}
    addParticles();

}

void Mantle::updateParams(void) {

    // override defaults with user parameters
    r->G = reb_sim_params.G;

    switch (reb_sim_params.config.gravity) {
        case gravity_solvers::TREE:
            r->gravity = reb_simulation::REB_GRAVITY_TREE;
            r->opening_angle2 = reb_sim_params.config.opening_angle2;
            r->softening = reb_sim_params.config.softening;
            break;
        case gravity_solvers::COMPENSATED:
            r->gravity = reb_simulation::REB_GRAVITY_COMPENSATED;
            break;
        case gravity_solvers::JACOBI:
            r->gravity = reb_simulation::REB_GRAVITY_JACOBI;
            break;
        case gravity_solvers::NONE:
            r->gravity = reb_simulation::REB_GRAVITY_NONE;
            break;
        default:
            break;
    }

    switch (reb_sim_params.config.collision) {
        case collision_solvers::DIRECT:
            r->collision = reb_simulation::REB_COLLISION_DIRECT;
            r->collision_resolve = reb_collision_resolve_hardsphere;
            break;
        case collision_solvers::LINE:
            r->collision = reb_simulation::REB_COLLISION_LINE;
            r->collision_resolve = reb_collision_resolve_hardsphere;
            break;
        case collision_solvers::LINETREE:
            r->collision = reb_simulation::REB_COLLISION_LINETREE;
            r->collision_resolve = reb_collision_resolve_hardsphere;
            break;
        case collision_solvers::TREE:
            r->collision = reb_simulation::REB_COLLISION_TREE;
            r->collision_resolve = reb_collision_resolve_hardsphere;
            break;
        default:
            break;
    }

    switch (reb_sim_params.config.integrator) {
        case integrators::WHFAST:
            r->integrator = reb_simulation::REB_INTEGRATOR_WHFAST;
            break;
        case integrators::SEI:
            r->integrator = reb_simulation::REB_INTEGRATOR_SEI;
            break;
        case integrators::LEAPFROG:
            r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG;
            break;
        default:
            break;
    }

    if ((reb_sim_params.config.gravity == gravity_solvers::TREE)
            || (reb_sim_params.config.collision == collision_solvers::TREE)
            || (reb_sim_params.config.collision == collision_solvers::LINETREE)) {
        reb_simulation_configure_box(r, reb_sim_params.config.box_size,
                reb_sim_params.config.x_boxes, reb_sim_params.config.y_boxes,
                reb_sim_params.config.z_boxes);
#ifdef OPENMP
        omp_set_num_threads(reb_sim_params.config.threads);
#endif
    }

    if (reb_sim_params.config.dt > 0) {
        r->dt = reb_sim_params.config.dt;
    }

    if (reb_sim_params.config.active > 0) {
        r->N_active = reb_sim_params.config.active;
    }

}

struct VisitState {

    void operator()(m_cart_state_vector state) {

        struct reb_particle p = {0};

        p.x = state.x; p.y = state.y; p.z = state.z;
        p.vx = state.vx; p.vy = state.vy; p.vz = state.vz;
        p.m = m;
        p.r = rp;
        p.hash = meta.hash;
        p.ap = const_cast<void*>(static_cast<const void*>(&meta));

        reb_simulation_add(r, p);

    }

    void operator()(m_kep_state_vector state) {

        struct reb_particle p = reb_particle_from_orbit(r->G, r->particles[0], m, state.a, state.e, state.inc, state.omega, state.Omega, state.f);

        p.r = rp;
        p.hash = meta.hash;
        p.ap = const_cast<void*>(static_cast<const void*>(&meta));

        reb_simulation_add(r, p);

    }

    reb_simulation* r;
    double m;
    double rp;
    const Meta& meta;

};

void Mantle::addParticles() {

    double LU = reb_sim_params.LU;
    double TU = reb_sim_params.TU;
    double MU = reb_sim_params.MU;

    std::partial_sort(particles.begin(), particles.begin()+1, particles.end(), [](auto p1, auto p2) {
            return p1.getMass() > p2.getMass();
            });

    for (Particle &p : particles) {

        p.setChar(LU, TU, MU);
        p.nonDim();
        std::visit(VisitState{r, p.getMass(), p.getRadius(), p.getMeta()}, p.getState());

    }

}

void Mantle::addHash() {

    uint32_t tmp_hash = 0;
    for (Particle &p : particles) {

        p.setHash(tmp_hash);

        tmp_hash++;

    }

}

void Mantle::addAdditionalForce() {

    r->additional_forces = m_additionalForces;
    r->force_is_velocity_dependent = 1;

}

void Mantle::addCollision(m_coll_context a) {

    r->collision_resolve = m_collisionResolve;
    std::for_each(particles.begin(), particles.end(), [=](auto& p) {p.setCollisionContext(a);});

}

void Mantle::addOde() {

    std::for_each(particles.begin(), particles.end(), [r=this->r] (Particle& p) {
            liftOptional(p.getOdes(), [r] (auto odes) {
                std::for_each(odes.get().begin(), odes.get().end(), [r] (const Ode& ode) {
                        struct reb_ode* ho = reb_ode_create(r, ode.y.size());
                        ho->derivatives = m_ode;
                        ho->ref = const_cast<void*>(static_cast<const void*>(&ode));
                        });
                });
            });

}

m_history Mantle::evalSim(double tf, int nt) {

    double dt = tf / (nt * reb_sim_params.TU);
    double t = 0.0;

    std::vector<std::vector<ParticleHistory>> state_history;
    for (int i = 0; i <= nt; i++) {
        reb_simulation_integrate(r, t);
        std::vector<ParticleHistory> tmp_ps(r->N);
        std::transform(r->particles, r->particles + r->N, tmp_ps.begin(), [t_=r->t] (auto p) {ParticleHistory result(p, t_); result.dim(); return result;});

        int count = std::count_if(tmp_ps.begin(), tmp_ps.end(), [](auto p) {return p.collectHistory();});
        std::vector<ParticleHistory> ps(count);
        std::copy_if(tmp_ps.begin(), tmp_ps.end(), ps.begin(), [](auto p) {return p.collectHistory();});
        state_history.push_back(ps);
        t += dt;
    }

    return state_history;

}

void m_additionalForces(struct reb_simulation* const r) {

    std::for_each(r->particles, r->particles + r->N, [t=r->t] (struct reb_particle& rp) {
            ParticleRestricted p(rp, t);
            liftOptional(p.getAdditionalForces(), [&] (std::vector<m_force_func>& as) {
                p.dim();
                auto accel = std::accumulate(as.begin(), as.end(), m_accel(0), [&](auto a, auto f) {return a + f(ParticleLite(p));});
                accel *= p.meta_props.TU * p.meta_props.TU / p.meta_props.LU;

                p.reb_p->ax += accel.ax;
                p.reb_p->ay += accel.ay;
                p.reb_p->az += accel.az;

            });

        });

}

int m_collisionResolve(struct reb_simulation* r, struct reb_collision c) {

    ParticleRestricted p1(r->particles[c.p1], r->t);
    ParticleRestricted p2(r->particles[c.p2], r->t);

    p1.dim();
    p2.dim();

    ParticleLite p1l(p1);
    ParticleLite p2l(p2);

    auto g = p1.getCollisionContext();
    auto f = g(p1l, p2l);
    int indx = f(p1l, p2l);

    p1.nonDim();
    p2.nonDim();

    p1.updateRebParticle();
    p2.updateRebParticle();

    p1.incCollision();
    p2.incCollision();
    r->collisions_log_n++;

    return indx;

}

void m_ode(struct reb_ode* const ode, double* const yDot, const double* const y, const double t) {

    auto pair = *static_cast<Ode*>(ode->ref);

    uint32_t hash = pair.hash;
    auto f = pair.ode_func;

    reb_particle* s = reb_simulation_particle_by_hash(ode->r, hash);

    if (s == NULL) {return;};

    ParticleRestricted p(*s, t);
    p.dim();

    ParticleLite pl(p);

    f(yDot, pl);

}

