#include "mantle_particle_restricted.h"
extern "C" {
	#include "rebound.h"
}

ParticleRestricted::ParticleRestricted(struct reb_particle &p, double t_) : meta_props(*static_cast<Meta*>(p.ap)) {

    state = {.x = p.x, .y = p.y, .z = p.z, .vx = p.vx, .vy = p.vy, .vz = p.vz, .t = t_};

    reb_p = &p;

}

const std::any& ParticleRestricted::getParticleProperties() const {return meta_props.particle_properties;}

const m_coll_context& ParticleRestricted::getCollisionContext() const {return meta_props.collisionContext;}

const std::optional<std::vector<m_force_func>>& ParticleRestricted::getAdditionalForces() const {return meta_props.additional_forces;}

void ParticleRestricted::updateRebParticle() {

    if (!stateUpdated) {return;}

    reb_p->x = state.x;
    reb_p->y = state.y;
    reb_p->z = state.z;
    reb_p->vx = state.vx;
    reb_p->vy = state.vy;
    reb_p->vz = state.vz;

}

void ParticleRestricted::nonDim() {

    Scale<std::divides<double>> dim_((std::divides<double>()), meta_props.LU, meta_props.TU);
    dim_(state);

}

void ParticleRestricted::dim() {

    Scale<std::multiplies<double>> dim_((std::multiplies<double>()), meta_props.LU, meta_props.TU);
    dim_(state);

}

ParticleLite::ParticleLite(ParticleRestricted &p_) : p(p_), state(p_.state) {

}

void ParticleLite::setState(double x, double y, double z, double vx, double vy, double vz) {

    p.stateUpdated = true;

    state.x = x;
    state.y = y;
    state.z = z;
    state.vx = vx;
    state.vy = vy;
    state.vz = vz;

}

ParticleHistory::ParticleHistory() {

}

ParticleHistory::ParticleHistory(struct reb_particle &p, double t_) {

    state = {.x = p.x, .y = p.y, .z = p.z, .vx = p.vx, .vy = p.vy, .vz = p.vz, .t = t_};

    meta_props = *static_cast<Meta*>(p.ap);

}

void ParticleHistory::nonDim() {

    Scale<std::divides<double>> dim_((std::divides<double>()), meta_props.LU, meta_props.TU);
    dim_(state);

}

void ParticleHistory::dim() {

    Scale<std::multiplies<double>> dim_((std::multiplies<double>()), meta_props.LU, meta_props.TU);
    dim_(state);

}

