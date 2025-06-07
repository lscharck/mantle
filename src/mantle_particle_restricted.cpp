#include "mantle_particle_restricted.hpp"

void ParticleRestricted::nonDim() {

    Scale<std::divides<double>> dim_((std::divides<double>()), meta_props.LU, meta_props.TU);
    dim_(state);

}

void ParticleRestricted::dim() {

    Scale<std::multiplies<double>> dim_((std::multiplies<double>()), meta_props.LU, meta_props.TU);
    dim_(state);

}

void ParticleLite::setState(double x, double y, double z, double vx, double vy, double vz) {

    state.x = x;
    state.y = y;
    state.z = z;
    state.vx = vx;
    state.vy = vy;
    state.vz = vz;

}

ParticleHistory::ParticleHistory(m_cart_state_vector _state, Meta m) {

    state = _state;

    meta_props = m;
    meta_props.syncOde();

}

void ParticleHistory::nonDim() {

    Scale<std::divides<double>> dim_((std::divides<double>()), meta_props.LU, meta_props.TU);
    dim_(state);

}

void ParticleHistory::dim() {

    Scale<std::multiplies<double>> dim_((std::multiplies<double>()), meta_props.LU, meta_props.TU);
    dim_(state);

}

