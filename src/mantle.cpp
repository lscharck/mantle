#include "mantle.hpp"
extern "C" {
	#include "rebound.h"
}

template<typename T, typename F>
void liftOptional(std::optional<T> val, F f) {
    if (val.has_value()) {f(val.value());}
}

Mantle::Mantle(const std::vector<std::unique_ptr<struct MConfig>>& user_configs, std::vector<Particle> userParticles) {

    r = reb_simulation_create();

    for (const auto& config: user_configs) {
        config->updateSim(this, r);
    }

    particles = userParticles;

    addHash();
    if (collisionContext.has_value()) {addCollision(collisionContext.value());};
    if (has_additional_forces) {addAdditionalForce();}
    if (has_ode) {addOde();}
    addParticles();

}

Mantle::~Mantle() {

    reb_simulation_free(r);

}

long Mantle::getCollisionLog() const {

    return r->collisions_log_n;

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

#ifdef OPENMP
    r->additional_forces = m_additionalForcesPar;
#else
    r->additional_forces = m_additionalForces;
#endif
    r->force_is_velocity_dependent = 1;

}

void Mantle::addCollision(m_coll_context a) {

    r->collision_resolve = m_collisionResolve;
    std::for_each(particles.begin(), particles.end(), [=](auto& p) {p.setCollisionContext(a);});

}

void Mantle::addOde() {

    std::for_each(particles.begin(), particles.end(), [TU=TU, r=this->r] (Particle& p) {
            liftOptional(p.getOdes(), [TU=TU, r] (auto odes) {
                std::for_each(odes.get().begin(), odes.get().end(), [TU=TU, r] (Ode& ode) {
                        struct reb_ode* ho = reb_ode_create(r, ode.size);
                        ho->derivatives = m_ode;
                        ode.TU = TU;
                        ode.init_y = Ode::nonDim(ode.init_y, 1.0, TU);
                        for (uint32_t i=0; i<ode.size;i++) {ho->y[i] = ode.init_y[i];}
                        ode.setOdey(ho);
                        ho->ref = const_cast<void*>(static_cast<const void*>(&ode));
                        });
                });
            });

}

m_history Mantle::evalSim(double tf, int nt) {

    double dt = tf / (nt * TU);
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

void m_additionalForcesPar(struct reb_simulation* const r) {

#pragma omp parallel for
    for (uint32_t i = 0; i < r->N; i++) {
            ParticleRestricted p(r->particles[i], r->t);
            liftOptional(p.getAdditionalForces(), [&] (std::vector<m_force_func>& as) {
                p.dim();
                auto accel = std::accumulate(as.begin(), as.end(), m_accel(0), [&](auto a, auto f) {return a + f(ParticleLite(p));});
                auto scale_accel = accel * (p.meta_props.TU * p.meta_props.TU / p.meta_props.LU);

                p.reb_p->ax += scale_accel.ax;
                p.reb_p->ay += scale_accel.ay;
                p.reb_p->az += scale_accel.az;

            });

        }

}

void m_additionalForces(struct reb_simulation* const r) {

    std::for_each(r->particles, r->particles + r->N, [t=r->t] (struct reb_particle& rp) {
            ParticleRestricted p(rp, t);
            liftOptional(p.getAdditionalForces(), [&] (std::vector<m_force_func>& as) {
                p.dim();
                auto accel = std::accumulate(as.begin(), as.end(), m_accel(0), [&](auto a, auto f) {return a + f(ParticleLite(p));});
                auto scale_accel = accel * (p.meta_props.TU * p.meta_props.TU / p.meta_props.LU);

                p.reb_p->ax += scale_accel.ax;
                p.reb_p->ay += scale_accel.ay;
                p.reb_p->az += scale_accel.az;

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

    auto _ode = *static_cast<Ode*>(ode->ref);

    uint32_t hash = _ode.hash;
    auto f = _ode.ode_func;

    reb_particle* s = reb_simulation_particle_by_hash(ode->r, hash);

    std::vector<double> _y;
    std::vector<double> y_dot;
    int size = _ode.size;

    if (s == NULL) {return;};

    ParticleRestricted p(*s, t);
    p.dim();

    ParticleLite pl(p);

    double TU = _ode.TU;
    _y.insert(_y.end(), y, y+size);
    y_dot = Ode::nonDim(f(Ode::dim(_y, 1.0, TU), pl), TU, TU*TU);
    for (int i = 0; i < size; i++) {
        yDot[i] = y_dot[i];
    }

}

