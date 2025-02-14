#include <iostream>
#include <fstream>
#include <functional>
#include <cassert>
#include <cmath>
#include "mantle.h"

struct properties1 {
    double param1;
    double param2;
    double mass;
    double radius;
};

struct properties2 {
    double param1;
    int param2;
    double mass;
    double radius;
};

m_accel f11(const ParticleLite& p1) {
    auto props_1 = p1.getParticleProperties<properties1>();
    double param1 = props_1->param1;
    double param2 = props_1->param2;
    double mass = props_1->mass;
    double radius = props_1->radius;
    assert((param1 == 17.8));
    assert((param2 == 90.3));
    assert((mass == 1.0));
    assert((radius == 0.4));

    return m_accel(0.01, 0.02, 0.03);
}

m_accel f12(const ParticleLite& p1) {
    auto props_1 = p1.getParticleProperties<properties1>();
    double param1 = props_1->param1;
    double param2 = props_1->param2;
    double mass = props_1->mass;
    double radius = props_1->radius;
    assert((param1 == 17.8));
    assert((param2 == 90.3));
    assert((mass == 1.0));
    assert((radius == 0.4));
    return m_accel(0.01, 0.01, 0.01);
}

m_accel f21(const ParticleLite& p2) {
    auto props_2 = p2.getParticleProperties<properties2>();
    double param1 = props_2->param1;
    int param2 = props_2->param2;
    double mass = props_2->mass;
    double radius = props_2->radius;
    assert((param1 = 45.2));
    assert((param2 = 1));
    assert((mass == 2.0));
    assert((radius == 0.8));
    return m_accel(0.01, 0.01, 0.01);
}

m_accel f22(const ParticleLite& p2) {
    auto props_2 = p2.getParticleProperties<properties2>();
    double param1 = props_2->param1;
    int param2 = props_2->param2;
    double mass = props_2->mass;
    double radius = props_2->radius;
    assert((param1 = 45.2));
    assert((param2 = 1));
    assert((mass == 2.0));
    assert((radius == 0.8));
    return m_accel(0.01, 0.01, 0.01);
}

void g11(double* const yDot,const ParticleLite& p1) {
    auto props_1 = p1.getParticleProperties<properties1>();
    double param1 = props_1->param1;
    double param2 = props_1->param2;
    double mass = props_1->mass;
    double radius = props_1->radius;
    assert((param1 == 17.8));
    assert((param2 == 90.3));
    assert((mass == 1.0));
    assert((radius == 0.4));
    yDot[0] = 0.1;
    yDot[1] = 0.2;
}

void g12(double* const yDot,const ParticleLite& p1) {
    auto props_1 = p1.getParticleProperties<properties1>();
    double param1 = props_1->param1;
    double param2 = props_1->param2;
    double mass = props_1->mass;
    double radius = props_1->radius;
    assert((param1 == 17.8));
    assert((param2 == 90.3));
    assert((mass == 1.0));
    assert((radius == 0.4));
    yDot[0] = 0.1;
    yDot[1] = 0.2;
    yDot[3] = 0.3;
}

void g21(double* const yDot,const ParticleLite& p2) {
    auto props_2 = p2.getParticleProperties<properties2>();
    double param1 = props_2->param1;
    int param2 = props_2->param2;
    double mass = props_2->mass;
    double radius = props_2->radius;
    assert((param1 = 45.2));
    assert((param2 = 1));
    assert((mass == 2.0));
    assert((radius == 0.8));
    yDot[0] = 0.1;
    yDot[1] = 0.2;
}

void g22(double* const yDot,const ParticleLite& p2) {
    auto props_2 = p2.getParticleProperties<properties2>();
    double param1 = props_2->param1;
    int param2 = props_2->param2;
    double mass = props_2->mass;
    double radius = props_2->radius;
    assert((param1 = 45.2));
    assert((param2 = 1));
    assert((mass == 2.0));
    assert((radius == 0.8));
    yDot[0] = 0.1;
    yDot[1] = 0.2;
    yDot[3] = 0.3;
}

std::vector<Particle> buildParticles() {
    double mass1 = 1.0;
    double radius1 = 0.4;
    Particle p1(mass1, radius1);
    p1.setStateCart(10, 0, 0, -2, 0, 0);
    p1.setCollisionTracker();
    p1.setHistoryTracker();
    properties1 props1 = {.param1 = 17.8, .param2 = 90.3, .mass=mass1, .radius=radius1};
    p1.setParticleProperties(props1);
    p1.setAdditionalForces({f11, f12});
    p1.setODE(g11, {0.0, 0.0});
    p1.setODE(g12, {0.0, 0.0, 0.0});

    double mass2 = 2.0;
    double radius2 = 0.8;
    Particle p2(mass2, radius2);
    p2.setStateCart(-10, 0, 0, 2, 0, 0);
    p2.setCollisionTracker();
    p2.setHistoryTracker();
    properties2 props2 = {.param1 = 45.2, .param2 = 1, .mass=mass2, .radius=radius2};
    p2.setParticleProperties(props2);
    p2.setAdditionalForces({f21, f22});
    p2.setODE(g21, {0.0, 0.0});
    p2.setODE(g22, {0.0, 0.0, 0.0});

    std::vector<Particle> ps = {p1, p2};

    return ps;
}

int elastic1D(ParticleLite& p1, ParticleLite& p2) {

    auto pr1 = p1.getParticleProperties<properties2>();
    double m1 = pr1->mass;
    auto pr2 = p2.getParticleProperties<properties1>();
    double m2 = pr2->mass;

    double v1 = p1.state.vx;
    double v2 = p2.state.vx;


    double v21 = (m1-m2)/(m1+m2)*v1 + 2*m2/(m1+m2)*v2;
    double v22 = 2*m1/(m1+m2)*v1 + (m2-m1)/(m1+m2)*v2;

    double dx1 = pr1->radius / 10;
    double dx2 = pr2->radius / 10;

    p1.setState(p1.state.x-dx1, p1.state.y, p1.state.z, v21, p1.state.vy, p1.state.vz);
    p2.setState(p2.state.x+dx2, p2.state.y, p2.state.z, v22, p2.state.vy, p2.state.vz);

    return 0;

}

m_coll_func decideCollision(const ParticleLite& p1, const ParticleLite& p2) {

    return elastic1D;

}

void saveFile(m_history ps) {

    std::fstream out_file;

    out_file.open("data.txt", std::ios::out);

    for (auto pt: ps) {
        for (auto p: pt) {
            auto s = p.state;
            std::print(out_file, "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} ", s.t, s.x, s.y, s.z, s.vx, s.vy, s.vz);
        }
        std::print(out_file, "\n");
    }

    out_file.close();

}


int full_test_setup() {

    m_sim_params params = {0};
    params.LU = 10.0;
    params.MU = 3.0;
    params.TU = 10.0;
    params.gravity = false;
    params.collision = true;

    Mantle m(params, buildParticles);

    m.setAdditionalForce();

    m.setCollisionContext(decideCollision);

    m.setOde();

    m.initMantle();

    auto hs = m.evalSim(10, 500);

    std::cout << m.getCollisionLog();
    saveFile(hs);

    return 0;

}

int main() {
    full_test_setup();
}

