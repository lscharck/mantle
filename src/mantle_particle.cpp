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
    ode.init_y   = y;
    ode.size     = y.size();
    ode.save_y   = ode.init_y;
    ode.id       = id;
    ode.particle_properties = particle_properties;

    if (!odes.has_value()) {
        odes = std::vector<Ode>{};
    }

    odes.value().push_back(ode);

}

void Meta::syncOde() {

    if (!odes.has_value()) {return;}

    for (auto& ode: odes.value()) {
        int size = ode.size;
        std::vector<double> tmp_y;
        tmp_y.insert(tmp_y.end(), ode.ho->y, ode.ho->y+size);
        ode.save_y = tmp_y;
    }

}

const std::optional<std::vector<double>> Meta::getOdey(uint32_t id) const {

    if (!odes.has_value()) { return std::nullopt;}

    for (auto ode: odes.value()) {
        if (ode.id == id) {
            return Ode::dim(ode.save_y, 1.0, TU);
        }

    }

    return std::nullopt;

}

template <typename Op>
std::vector<double> scale(Op op, std::vector<double> y, double s1, double s2) {

    std::vector<double> new_y = y;

    for (uint32_t i=0; i<y.size(); i+=2) {
        new_y[i] = op(y[i], s1);
    }

    for (uint32_t i=1; i<y.size(); i+=2) {
        new_y[i] = op(y[i], s2);
    }

    return new_y;

}

std::vector<double> Ode::dim(std::vector<double> y, double s1, double s2) {
    return scale((std::divides<double>()), y, s1, s2);
}

std::vector<double> Ode::nonDim(std::vector<double> y, double s1, double s2) {
    return scale((std::multiplies<double>()), y, s1, s2);
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

void Particle::nonDim() {

    Scale<std::divides<double>> dim_((std::divides<double>()), meta_props.LU, meta_props.TU);
    std::visit(dim_, state);

    m = m / meta_props.MU;
    r = r / meta_props.LU;

}

void Particle::dim() {

    Scale<std::multiplies<double>> dim_((std::multiplies<double>()), meta_props.LU, meta_props.TU);
    std::visit(dim_, state);

    m = m * meta_props.MU;
    r = r * meta_props.LU;

}

m_accel operator+(const m_accel& t, const m_accel& u) {
    return {.ax = t.ax + u.ax, .ay = t.ay + u.ay, .az = t.az + u.az};
}

m_accel operator*(double s, const m_accel& u) {
    return {.ax = s*u.ax, .ay = s*u.ay, .az = s*u.az};
}

m_accel operator*(const m_accel& u, double s) {
    return {.ax = s*u.ax, .ay = s*u.ay, .az = s*u.az};
}

