#include "mantle_particle.h"

Meta::Meta() {

    collisions = -1;
    count_collisions = false;
    collect_history = false;
    hash = 0;
    LU = 1;
    TU = 1;
    MU = 1;

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

void Particle::setCollisionTracker(void) {

    meta_props.count_collisions = true;

}

void Particle::setHistoryTracker(void) {

    meta_props.collect_history = true;

}

void Particle::setAdditionalForces(std::vector<m_force_func> function_array) {

    meta_props.additional_forces = function_array;

}

void Particle::setODE(m_ode_func func, std::vector<double> y) {

    Ode ode;
    ode.ode_func = func;
    ode.y = y;
    ode.particle_properties = meta_props.particle_properties;
    ode.id = 0;

    if (!meta_props.odes.empty()) {ode.id = meta_props.odes.back().id++;}

    meta_props.odes.push_back(ode);

}

void Particle::setCollisionContext(m_coll_context collision_context) {meta_props.collisionContext = collision_context;}

void Particle::setChar(double LU, double TU, double MU) {meta_props.LU = LU; meta_props.TU = TU; meta_props.MU = MU;};

void Particle::setHash(uint32_t hash) {meta_props.hash = hash;};

const std::variant<m_cart_state_vector, m_kep_state_vector>& Particle::getState() const {return state;}

const Meta& Particle::getMeta() const {return meta_props;}

const std::vector<Ode>& Particle::getODE() const {return meta_props.odes;}

double Particle::getMass() const {return m;}

double Particle::getRadius() const {return r;}

uint32_t Particle::getHash() const {return meta_props.hash;};

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

m_accel operator*=(double s, const m_accel& u) {
    return {.ax = s*u.ax, .ay = s*u.ay, .az = s*u.az};
}

m_accel operator*=(const m_accel& u, double s) {
    return {.ax = s*u.ax, .ay = s*u.ay, .az = s*u.az};
}

