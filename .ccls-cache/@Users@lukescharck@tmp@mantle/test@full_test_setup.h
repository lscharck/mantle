#include <iostream>
#include <functional>
#include <cassert>
#include <mantle.h>

struct properties1 {
    double param1;
    double param2;
};

struct properties2 {
    double param1;
    int param2;
};

m_accel f11(ParticleLite p) {
    double x = p.x;
    double y = p.y;
    double z = p.z;
    double vx = p.vx;
    double vy = p.vy;
    double vz = p.vz;

    properties1* props = (properties1*)p.getParticleProperties();
    double pr1 = props->param1;
    double pr2 = props->param2;
    assert(pr1 == 17.8);
    assert(pr2 == 90.3);

    return m_accel(0.01, 0.02, 0.03);
}

m_accel f12(ParticleLite p) {
    return m_accel(0.01, 0.01, 0.01);
}

m_accel f21(ParticleLite p) {
    properties2* props = (properties2*)p.getParticleProperties();
    double pr1 = props->param1;
    int pr2 = props->param2;
    assert(pr1 = 45.2);
    assert(pr2 = 1);
    return m_accel(0.01, 0.01, 0.01);
}

m_accel f22(ParticleLite p) {
    return m_accel(0.01, 0.01, 0.01);
}

void g11(const double* const y, double* const yDot, const double t, void* p) {
    properties1* pr = (properties1*)p;
    double pr1 = pr->param1;
    double pr2 = pr->param2;
    assert(pr1 == 17.8);
    assert(pr2 == 90.3);
    yDot[0] = 0.1;
    yDot[1] = 0.2;
}

void g12(const double* const y, double* const yDot, const double t, void* p) {
    properties1* pr = (properties1*)p;
    double pr1 = pr->param1;
    double pr2 = pr->param2;
    assert(pr1 == 17.8);
    assert(pr2 == 90.3);
    yDot[0] = 0.1;
    yDot[1] = 0.2;
    yDot[3] = 0.3;
}

void g21(const double* const y, double* const yDot, const double t, void* p) {
    properties2* pr = (properties2*)p;
    double pr1 = pr->param1;
    int pr2 = pr->param2;
    assert(pr1 = 45.2);
    assert(pr2 = 1);
    yDot[0] = 0.1;
    yDot[1] = 0.2;
}

void g22(const double* const y, double* const yDot, const double t, void* p) {
    properties2* pr = (properties2*)p;
    double pr1 = pr->param1;
    int pr2 = pr->param2;
    assert(pr1 = 45.2);
    assert(pr2 = 1);
    yDot[0] = 0.1;
    yDot[1] = 0.2;
    yDot[3] = 0.3;
}

std::vector<Particle> buildParticles() {
    double mass1 = 1.0;
    double radius1 = 0.5;
    Particle p1(mass1, radius1);
    p1.setStateCart(10, 0, 0, -2, 0, 0);
    p1.setCollisionTracker();
    p1.setHistoryTracker();
    properties1 props1 = {.param1 = 17.8, .param2 = 90.3};
    p1.setUserProperties(&props1, sizeof(properties1));
    p1.setAdditionalForces({f11, f12});
    p1.setODE(g11, {0.0, 0.0});
    p1.setODE(g12, {0.0, 0.0, 0.0});

    double mass2 = 2.0;
    double radius2 = 0.7;
    Particle p2(mass2, radius2);
    p2.setStateCart(-10, 0, 0, 2, 0, 0);
    p2.setCollisionTracker();
    p2.setHistoryTracker();
    properties2 props2 = {.param1 = 45.2, .param2 = 1};
    p2.setUserProperties(&props2, sizeof(properties2));
    p2.setAdditionalForces({f21, f22});
    p1.setODE(g21, {0.0, 0.0});
    p1.setODE(g22, {0.0, 0.0, 0.0});

    std::vector<Particle> ps = {p1, p2};

    return ps;
}

int f(ParticleLite& p1, ParticleLite& p2) {

    return -1;

}

m_coll_func decideCollision(ParticleLite p1, ParticleLite p2) {

    return f;

}


int full_test_setup() {

    m_sim_params params = {0};
    params.LU = 1;
    params.MU = 1;
    params.TU = 1;
    params.G = 1;

    Mantle m(params, buildParticles);

    m.setAdditionalForce();

    m.setCollisionContext(decideCollision);

    m.setOde();

    m.initMantle();

    auto hs = m.evalSim(100, 50);

}
