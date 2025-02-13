#include "mantle_particle_restricted.h"
extern "C" {
	#include "rebound.h"
}

ParticleRestricted::ParticleRestricted() {

}

ParticleRestricted::ParticleRestricted(struct reb_particle &p, double t_) {

    x = p.x;
    y = p.y;
    z = p.z;
    vx = p.vx;
    vy = p.vy;
    vz = p.vz;
    t = t_;

    meta_props = (Meta*)p.ap;
    reb_p = &p;
    LU = meta_props->LU;
    TU = meta_props->TU;
    VU = LU/TU;

}


void ParticleRestricted::setState(double x_, double y_, double z_, double vx_, double vy_, double vz_) {

    stateUpdated = true;

    x = x_;
    y = y_;
    z = z_;
    vx = vx_;
    vy = vy_;
    vz = vz_;

}

const std::any& ParticleRestricted::getParticleProperties() const {return meta_props->user_properties;}

const m_coll_context& ParticleRestricted::getCollisionContext() const {return meta_props->collisionContext;}

const std::optional<std::vector<m_force_func>>& ParticleRestricted::getAdditionalForces() const {return meta_props->additional_forces;}

void ParticleRestricted::updateRebParticle() {

    if (!stateUpdated) {return;}

    reb_p->x = x;
    reb_p->y = y;
    reb_p->z = z;
    reb_p->vx = vx;
    reb_p->vy = vy;
    reb_p->vz = vz;

}

template <typename Op>
void ParticleRestricted::scale(Op op) {

    x = op(x, LU);
    y = op(y, LU);
    z = op(z, LU);
    vx = op(vx, VU);
    vy = op(vy, VU);
    vz = op(vz, VU);
    t = op(t, TU);

}

void ParticleRestricted::nonDim() {scale((std::divides<double>()));}

void ParticleRestricted::dim() {scale((std::multiplies<double>()));}

ParticleLite::ParticleLite(ParticleRestricted &p_) {

    x = p_.x;
    y = p_.y;
    z = p_.z;
    vx = p_.vx;
    vy = p_.vy;
    vz = p_.vz;

    p = &p_;

}

void ParticleLite::setState(double x_, double y_, double z_, double vx_, double vy_, double vz_) {p->setState(x_, y_, z_, vx_, vy_, vz_);}

