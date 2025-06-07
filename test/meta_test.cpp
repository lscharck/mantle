#include <iostream>
#include <fstream>
#include <functional>
#include <cassert>
#include <cmath>
#include "mantle.hpp"

std::vector<double> dummy_ode(std::vector<double> y, const ParticleLite& p1) {
    return {0.0};
}

m_accel dummy_forcing(const ParticleLite& p1) {
    return m_accel(0.01, 0.02, 0.03);
}

int dummy_collision(ParticleLite& p1, ParticleLite& p2) {
    return 0;
}

m_coll_func dummy_decide(const ParticleLite& p1, const ParticleLite& p2) {
    return dummy_collision;
}

struct dummy_properties {
    double param1;
    double param2;
    double mass;
    double radius;
};

bool cmpr(struct dummy_properties p1, struct dummy_properties p2) {
    return ((p1.param1 = p2.param1) && (p1.param2 = p2.param2) && (p1.mass = p2.mass) && (p1.radius = p2.radius));
}

bool setOde_test() {
    Meta meta;
    std::vector<double> y = {1.0, 2.0};
    uint32_t id = 2;
    struct dummy_properties p_in = {.param1=2.2, .param2=3.8, .mass=9.3, .radius=100.56};
    meta.particle_properties = p_in;
    meta.setOde(dummy_ode, y, id);
    Ode ode = meta.odes.value()[0];
    struct dummy_properties p_out = std::any_cast<struct dummy_properties>(ode.particle_properties);
    return ((ode.id == id) && (ode.ode_func == dummy_ode) && (ode.init_y == y) && (cmpr(p_in, p_out)));
}

bool setCollisionTracker_test() {
    Meta meta;
    meta.setCollisionTracker();
    return (meta.count_collisions == true);
}

bool setHistoryTracker_test() {
    Meta meta;
    meta.setHistoryTracker();
    return (meta.collect_history == true);
}

bool setAdditionalForces_test() {
    Meta meta;
    std::vector<m_force_func> fs = {dummy_forcing};
    meta.setAdditionalForces(fs);
    return (meta.additional_forces == fs);
}

bool setCollisionContext_test() {
    Meta meta;
    meta.setCollisionContext(dummy_decide);
    return (meta.collision_context == dummy_decide);
}

bool setChar_test() {
    Meta meta;
    double LU = 10.0;
    double TU = 20.0;
    double MU = 30.0;
    meta.setChar(LU, TU, MU);
    return ((meta.LU == LU) && (meta.TU == TU) && (meta.MU == MU));
}

bool setHash_test() {
    Meta meta;
    std::vector<double> y = {1.0, 2.0};
    meta.setOde(dummy_ode, y, 2);
    uint32_t hash_in = 2;
    meta.setHash(hash_in);
    uint32_t hash_out = meta.hash;
    uint32_t hash_ode_out = meta.odes.value()[0].hash;
    return ((hash_out == hash_in) && (hash_ode_out = hash_in));
}

//bool saveOde_test() {
//void saveOde();
//}

bool setParticleProperties_test() {
    Meta meta;
    struct dummy_properties p_in = {.param1=2.2, .param2=3.8, .mass=9.3, .radius=100.56};
    meta.setParticleProperties(p_in);
    struct dummy_properties p_out = std::any_cast<struct dummy_properties>(meta.particle_properties);
    return cmpr(p_in, p_out);
}

bool incCollision_test() {
    Meta meta;
    meta.count_collisions = true;
    meta.collisions = 0;
    meta.incCollision();
    return (meta.collisions == 1);
}

bool getCollisionLog_test() {
    Meta meta;
    int c_in = 1;
    meta.collisions = c_in;
    int c_out = meta.getCollisionLog();
    return (c_out == c_in);
}

bool collectHistory_test() {
    Meta meta;
    meta.collectHistory();
    return (meta.collect_history == true);
}

/*
bool getOdes_test() {
const std::optional<std::reference_wrapper<std::vector<Ode>>> getOdes() {if (odes.has_value()) {return std::ref(odes.value());} else {return std::nullopt;}}
}

bool getOdey_test() {
const std::optional<std::vector<double>> getOdey(uint32_t id) const;
}
*/

bool getHash_test() {
    Meta meta;
    uint32_t hash_in = 34;
    meta.hash = hash_in;
    uint32_t hash_out = meta.getHash();
    return (hash_in == hash_out);
}

bool meta_test() {

    std::string log = "";
    bool t1 = setOde_test();
    if (!t1) {log+="setOde_test failed";}

    bool t2 = setCollisionTracker_test();
    if (!t2) {log+=" setCollisionTracker_test failed";}

    bool t3 = setHistoryTracker_test();
    if (!t3) {log+=" setHistoryTracker_test failed";}

    bool t4 = setAdditionalForces_test();
    if (!t4) {log+=" setAdditionalForces_test failed";}

    bool t5 = setCollisionContext_test();
    if (!t5) {log+=" setCollisionContext_test failed";}

    bool t6 = setChar_test();
    if (!t6) {log+=" setChar_test failed";}

    bool t7 = setHash_test();
    if (!t7) {log+=" setHash_test failed";}

    bool t8 = setParticleProperties_test();
    if (!t8) {log+=" setParticleProperties_test failed";}

    bool t9 = incCollision_test();
    if (!t9) {log+=" incCollision_test failed";}

    bool t10 = getCollisionLog_test();
    if (!t10) {log+=" getCollisionLog_test failed";}

    bool t11 = collectHistory_test();
    if (!t11) {log+=" collectHistory_test failed";}

    bool t12 = getHash_test();
    if (!t12) {log+=" getHash_test failed";}

    std::cout << log << std::endl;

    return (t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12);

}

