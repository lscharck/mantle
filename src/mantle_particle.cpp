#include "mantle_particle.hpp"
extern "C" {
	#include "rebound.h"
}

Meta::Meta() {

    collisions       = 0;
    count_collisions = false;
    collect_history  = false;
    removable        = true;
    hash             = 0;
    LU               = 1.0;
    TU               = 1.0;
    MU               = 1.0;

}

void Meta::setOde(m_ode_func func, std::vector<double> y, uint32_t id) {

    Ode ode;
    ode.ode_func = func;
    ode.size     = y.size();
    ode.id       = id;
    ode.particle_properties = particle_properties;
    ode.particle_properties_gen = particle_properties_gen;

    for (uint32_t i=0; i<ode.size; i++) {ode.init_y[i] = y[i];};
    ode.save_y   = ode.init_y;

    if (!odes.has_value()) {
        odes = std::vector<Ode>{};
    }

    odes.value().push_back(ode);

}

void Meta::syncOde() {

    if (!odes.has_value()) {return;}

    for (auto& ode: odes.value()) {
        int size = ode.size;
        Ode::dim(ode.ho->y, 1.0, TU, ode.size);
        std::copy(ode.ho->y, ode.ho->y+size, ode.save_y.begin());
        Ode::nonDim(ode.ho->y, 1.0, TU, ode.size);
    }

}

void Meta::nonDimId(uint32_t id) {

    if (!odes.has_value()) {return;}

    for (auto& ode: odes.value()) {
        if (ode.id == id) {Ode::nonDim(ode.ho->y, 1.0, TU, ode.size);}
    }

}


const std::optional<std::array<double, MAX_STATE>> Meta::getSaveOdey(uint32_t id) const {

    if (!odes.has_value()) {return std::nullopt;}

    for (auto& ode: odes.value()) {
        if (ode.id == id) {return ode.save_y;}
    }

    return std::nullopt;

}

const std::optional<double*> Meta::getOdey(uint32_t id) const {

    if (!odes.has_value()) {return std::nullopt;}

    for (auto& ode: odes.value()) {
        if (ode.id == id) {Ode::dim(ode.ho->y, 1.0, TU, ode.size); return ode.ho->y;}
    }

    return std::nullopt;

}

Particle::Particle(double mass, double radius) : meta_props() {

    m = mass;
    r = radius;

}

void Particle::setStateCart(double x, double y, double z, double vx, double vy, double vz) {

    m_cart_state_vector s = {0};
    s.x = x;
    s.y = y;
    s.z = z;
    s.vx = vx;
    s.vy = vy;
    s.vz = vz;

    state = s;

}

void Particle::setStateKep(double a, double e, double inc, double Omega, double omega, double f) {

    m_kep_state_vector s = {0};
    s.a = a;
    s.e = e;
    s.inc = inc;
    s.Omega = Omega;
    s.omega = omega;
    s.f = f;

    state = s;

}

void Meta::setHash(uint32_t _hash) {

    hash = _hash;

    if (!odes.has_value()) {return;}

    for (auto &o: odes.value()) {
        o.hash = hash;
    }

}

