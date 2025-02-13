#include "mantle.h"
extern "C" {
	#include "rebound.h"
}

template<typename T, typename F>
void liftOptional(std::optional<T> val, F f) {
    if (val.has_value()) {f(val.value());}
}

// this function runs mantle with a user defines set up
Mantle::Mantle(m_sim_params user_params, userParticles userFn) {

    reb_sim_params = user_params;
    userParticlesFn = userFn;
    has_additional_forces = false;
    has_ode = false;

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

void Mantle::initMantle(void) {

    r = reb_simulation_create();
    updateParams();
    particles = userParticlesFn();
    if (collisionContext.has_value()) {addCollision(collisionContext.value());};
    if (has_additional_forces) {addAdditionalForce();}
    if (has_ode) {addOde();}
    addParticles();

}

// this function updates mantle and reb with user parameters
void Mantle::updateParams(void) {

    // override defaults with user parameters
    r->G = reb_sim_params.G;
    if (!reb_sim_params.gravity) {r->gravity = reb_simulation::REB_GRAVITY_NONE;};
    if (reb_sim_params.default_collision) {r->collision = reb_simulation::REB_COLLISION_LINE; r->collision_resolve = reb_collision_resolve_hardsphere;};
    if (reb_sim_params.collision) {r->collision = reb_simulation::REB_COLLISION_LINE;};

}

// this function preforms the addition of particle to reb
// adds meta parameters to meta
// adds additional forces and ode functions inside meta
struct VisitState {

    void operator()(m_cart_state_vector state) {

        struct reb_particle p = {0};

        p.x = state.x; p.y = state.y; p.z = state.z;
        p.vx = state.vx; p.vy = state.vy; p.vz = state.vz;
        p.m = m;
        p.r = rp;
        p.hash = meta->hash;
        p.ap = (void*)meta;

        reb_simulation_add(r, p);

    }

    void operator()(m_kep_state_vector state) {

        struct reb_particle p = reb_particle_from_orbit(r->G, r->particles[0], m, state.a, state.e, state.inc, state.omega, state.Omega, state.f);

        p.r = rp;
        p.hash = meta->hash;
        p.ap = (void*)meta;

        reb_simulation_add(r, p);

    }

    reb_simulation* r;
    double m;
    double rp;
    Meta* meta;

};

// this function add user defines particles to reb via a user callback
void Mantle::addParticles(void) {

    double LU = reb_sim_params.LU;
    double TU = reb_sim_params.TU;
    double MU = reb_sim_params.MU;

    std::partial_sort(particles.begin(), particles.begin()+1, particles.end(), [](auto p1, auto p2) {
            return p1.getMass() > p2.getMass();
            });

    uint32_t tmp_hash = 0;
    for (Particle &p : particles) {

        p.setChar(LU, TU, MU);
        p.nonDim();
        p.setHash(tmp_hash);
        std::visit(VisitState{r, p.getMass(), p.getRadius(), p.getMeta()}, p.getState());

        tmp_hash++;

    }

}

void Mantle::addAdditionalForce(void) {

    r->additional_forces = m_additionalForces;

}

void Mantle::addCollision(m_coll_context a) {

    r->collision_resolve = m_collisionResolve;
    std::for_each(particles.begin(), particles.end(), [=](auto& p) {p.setCollisionContext(a);});

}

// this function adds callbacks to reb
void Mantle::addOde(void) {

    std::for_each(particles.begin(), particles.end(), [r=this->r, &otp=this->ode_to_particle] (Particle& p) {
            std::for_each(p.getODE()->begin(), p.getODE()->end(), [&p, &otp, r] (Ode& ode) {
                    struct reb_ode* ho = reb_ode_create(r,ode.y.size());
                    ho->derivatives = m_ode;
                    //ho->ref = (void*)&ode;
                    otp.push_back(std::make_pair(ode.id, p.meta_props.hash));
                    ho->ref = (void*)&otp.back();
                    });
            });

}

// this function evaluates the sim collecting time history
std::vector<std::vector<ParticleRestricted>> Mantle::evalSim(double tf, int nt) {

    double dt = tf / nt;
    double t = dt;

    std::vector<std::vector<ParticleRestricted>> state_history;
    for (int i = 0; i < nt; i++) {
        reb_simulation_integrate(r, t);
        std::vector<ParticleRestricted> ps(r->N);
        std::transform(r->particles, r->particles + r->N, ps.begin(), [t_=r->t] (auto p) {ParticleRestricted result(p, t_); result.dim(); return result;});
        t += dt;
        state_history.push_back(ps);
    }

    return state_history;

}

void m_additionalForces(struct reb_simulation* const r) {

    std::for_each(r->particles, r->particles + r->N, [t=r->t] (auto rp) {
            ParticleRestricted p(rp, t);
            liftOptional(p.getAdditionalForces(), [&] (auto as) {
                p.dim();
                auto accel = std::accumulate(as.begin(), as.end(), m_accel(0), [&](auto a, auto f) {return a + f(ParticleLite(p));});
                accel *= p.TU * p.TU / p.LU;

                p.reb_p->ax = accel.ax;
                p.reb_p->ay = accel.ay;
                p.reb_p->az = accel.az;

            });

        });

}

// this function calls the specific collision function per particle
int m_collisionResolve(struct reb_simulation* r, struct reb_collision c) {

    struct reb_particle* s1 = &r->particles[c.p1];
    struct reb_particle* s2 = &r->particles[c.p2];
    ParticleRestricted p1(r->particles[c.p1], r->t);
    ParticleRestricted p2(r->particles[c.p2], r->t);

    p1.dim();
    p2.dim();

    ParticleLite p1l(p1);
    ParticleLite p2l(p2);

    m_coll_context g = p1.getCollisionContext();
    m_coll_func f = g(p1l, p2l);
    int indx = f(p1l, p2l);

    p1.nonDim();
    p2.nonDim();

    s1->x = p1.x;
    s1->y = p1.y;
    s1->z = p1.z;
    s1->vx = p1.vx;
    s1->vy = p1.vy;
    s1->vz = p1.vz;

    s2->x = p2.x;
    s2->y = p2.y;
    s2->z = p2.z;
    s2->vx = p2.vx;
    s2->vy = p2.vy;
    s2->vz = p2.vz;

    return indx;

}

// this function calls the specific ode function per particle
void m_ode(struct reb_ode* const ode, double* const yDot, const double* const y, const double t) {

    std::pair<uint32_t, uint32_t>* pair = (std::pair<uint32_t, uint32_t>*)ode->ref;

    uint32_t ode_id = std::get<0>(*pair);
    uint32_t hash = std::get<1>(*pair);

    reb_particle* s = reb_simulation_particle_by_hash(ode->r, hash);

    if (s == NULL) {return;};

    ParticleRestricted p(*s, t);
    p.dim();

    auto f = p.meta_props->odes[ode_id].ode_func;

    ParticleLite pl(p);

    f(yDot, pl);

}

