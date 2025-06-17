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
    r->extras = const_cast<void*>(static_cast<const void*>(this));

    for (const auto& config: user_configs) {
        config->updateSim(this, r);
    }

    particles = userParticles;

    addHash();
    if (collisionContext.has_value()) {addCollision(collisionContext.value());};
    if (has_additional_forces) {addAdditionalForce();}
    if (has_ode) {addOde();}
    if (removeParticleTest.has_value()) {addHeartbeat();}
    addParticles();

}

Mantle::~Mantle() {

    std::for_each(particles.begin(), particles.end(), [] (Particle& p) {
        liftOptional(p.getOdes(), [] (auto odes) {std::for_each(odes.get().begin(), odes.get().end(), [] (Ode& ode) {if (ode.ho != NULL) {reb_ode_free(ode.ho);}});});
        });
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
    if (force_is_velocity_dep) {r->force_is_velocity_dependent = 1;}

}

void Mantle::addHeartbeat() {

#ifdef OPENMP
    r->heartbeat = m_heartbeatPar;
#else
    r->heartbeat = m_heartbeat;
#endif

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
                        Ode::nonDim(ode.init_y.data(), 1.0, TU, ode.size);
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
    state_history.reserve(nt);
    for (int i = 0; i <= nt; i++) {
        reb_simulation_integrate(r, t);
        std::vector<ParticleHistory> tmp_ps(r->N);
        std::transform(r->particles, r->particles + r->N, tmp_ps.begin(), [t_=r->t] (auto p) {ParticleHistory result(m_rebToCart(p, t_), *static_cast<Meta*>(p.ap)); result.dim(); return result;});

        int count = std::count_if(tmp_ps.begin(), tmp_ps.end(), [](auto p) {return p.collectHistory();});
        std::vector<ParticleHistory> ps(count);
        std::copy_if(tmp_ps.begin(), tmp_ps.end(), ps.begin(), [](auto p) {return p.collectHistory();});
        state_history.emplace_back(ps);
        if (r->status == REB_STATUS::REB_STATUS_USER) {break;}
        t += dt;
    }

#ifdef MPI
    reb_mpi_finalize(r);
#endif

    return state_history;

}

void m_additionalForcesPar(struct reb_simulation* const r) {

#pragma omp parallel for
    for (uint32_t i = 0; i < r->N; i++) {
            ParticleRestricted p(m_rebToCart(r->particles+i, r->t), *static_cast<Meta*>(r->particles[i].ap));
            liftOptional(p.getAdditionalForces(), [&] (std::vector<m_force_func>& as) {
                p.dim();
                auto accel = std::accumulate(as.begin(), as.end(), m_accel(0), [&](auto a, auto f) {return a + f(ParticleLite(p));});
                auto scale_accel = accel * (p.meta_props.TU * p.meta_props.TU / p.meta_props.LU);

                r->particles[i].ax += scale_accel.ax;
                r->particles[i].ay += scale_accel.ay;
                r->particles[i].az += scale_accel.az;

            });

        }

}

void m_additionalForces(struct reb_simulation* const r) {

    std::for_each(r->particles, r->particles + r->N, [t=r->t] (struct reb_particle& rp) {
            ParticleRestricted p(m_rebToCart(rp, t), *static_cast<Meta*>(rp.ap));
            liftOptional(p.getAdditionalForces(), [&] (std::vector<m_force_func>& as) {
                p.dim();
                auto accel = std::accumulate(as.begin(), as.end(), m_accel(0), [&](auto a, auto f) {return a + f(ParticleLite(p));});
                auto scale_accel = accel * (p.meta_props.TU * p.meta_props.TU / p.meta_props.LU);

                rp.ax += scale_accel.ax;
                rp.ay += scale_accel.ay;
                rp.az += scale_accel.az;

            });

        });

}

int m_collisionResolve(struct reb_simulation* r, struct reb_collision c) {

    struct reb_particle* reb_p1 = r->particles+c.p1;
    struct reb_particle* reb_p2 = r->particles+c.p2;

    ParticleRestricted p1(m_rebToCart(reb_p1, r->t), *static_cast<Meta*>(reb_p1->ap));
    ParticleRestricted p2(m_rebToCart(reb_p2, r->t), *static_cast<Meta*>(reb_p2->ap));

    p1.dim();
    p2.dim();
    ParticleLite p1l(p1);
    ParticleLite p2l(p2);

    auto g = p1.getCollisionContext();
    auto f = g(p1l, p2l);
    int indx = f(p1l, p2l);

    p1.nonDim();
    p2.nonDim();

    auto s1 = p1.updateRebParticle();
    reb_p1->x = s1.x;
    reb_p1->y = s1.y;
    reb_p1->z = s1.z;
    reb_p1->vx = s1.vx;
    reb_p1->vy = s1.vy;
    reb_p1->vz = s1.vz;
    auto s2 = p2.updateRebParticle();
    reb_p2->x = s2.x;
    reb_p2->y = s2.y;
    reb_p2->z = s2.z;
    reb_p2->vx = s2.vx;
    reb_p2->vy = s2.vy;
    reb_p2->vz = s2.vz;

    p1.incCollision();
    p2.incCollision();
    r->collisions_log_n++;

    return indx;

}

void m_ode(struct reb_ode* const ode, double* const yDot, const double* const y, const double t) {

    auto _ode = static_cast<Ode*>(ode->ref);
    double TU = _ode->TU;
    uint32_t hash = _ode->hash;
    auto f = _ode->ode_func;

    reb_particle* s = reb_simulation_particle_by_hash(ode->r, hash);

    if (s == NULL) {return;};

    ParticleRestricted p(m_rebToCart(s, t), *static_cast<Meta*>(s->ap));
    p.dim();

    Ode::dim(const_cast<double*>(y), 1.0, TU, _ode->size);
    f(yDot, y, ParticleLite(p));

    Ode::nonDim(yDot, TU, TU*TU, _ode->size);
    Ode::nonDim(const_cast<double*>(y), 1.0, TU, _ode->size);

}

void m_heartbeatPar(struct reb_simulation* r) {

    auto m = static_cast<Mantle*>(r->extras);
    auto tester = m->removeParticleTest.value();
    auto terminate = m->simTerminate.value();

    std::vector<ParticleRestricted> particles_to_remove;

    #pragma omp parallel 
    {
        std::vector<ParticleRestricted> thread_local_ids;

        #pragma omp for nowait
        for (uint32_t i = 0; i < r->N; i++) {
            ParticleRestricted p(m_rebToCart(r->particles+i, r->t), *static_cast<Meta*>(r->particles[i].ap));
            p.dim();
            if (p.isRemovable() && tester(ParticleLite(p))) {
                thread_local_ids.push_back(p);
            }

        }

        #pragma omp critical
        {
            particles_to_remove.insert(particles_to_remove.end(), thread_local_ids.begin(), thread_local_ids.end());
        }

    }

    for (auto p : particles_to_remove) {
        liftOptional(p.getOdes(), [] (auto odes) {std::for_each(odes.get().begin(), odes.get().end(), [] (Ode& ode) {reb_ode_free(ode.ho);});});
        reb_simulation_remove_particle_by_hash(r, p.getHash(), 0);
    }

    bool yes_term = terminate(r->N);
    if (yes_term) {reb_simulation_stop(r);}

}

void m_heartbeat(struct reb_simulation* r) {

    auto m = static_cast<Mantle*>(r->extras);
    auto tester = m->removeParticleTest.value();
    auto terminate = m->simTerminate.value();

    std::for_each(r->particles, r->particles + r->N, [tester=tester, t=r->t, r=r] (struct reb_particle& rp) {
            ParticleRestricted p(m_rebToCart(rp, t), *static_cast<Meta*>(rp.ap));
            p.dim();
            if (p.isRemovable() && tester(ParticleLite(p))) {
                liftOptional(p.getOdes(), [] (auto odes) {std::for_each(odes.get().begin(), odes.get().end(), [] (Ode& ode) {reb_ode_free(ode.ho); ode.ho = NULL;});});
                reb_simulation_remove_particle_by_hash(r, p.getHash(), 0);
            }
    });

    bool yes_term = terminate(r->N);
    if (yes_term) {reb_simulation_stop(r);}

}

m_cart_state_vector m_rebToCart(struct reb_particle& p, double t_) {

    return {.x = p.x, .y = p.y, .z = p.z, .vx = p.vx, .vy = p.vy, .vz = p.vz, .t = t_};

}

m_cart_state_vector m_rebToCart(struct reb_particle* p, double t_) {

    return {.x = p->x, .y = p->y, .z = p->z, .vx = p->vx, .vy = p->vy, .vz = p->vz, .t = t_};

}

