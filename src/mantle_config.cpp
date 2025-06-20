#include "mantle_config.hpp"
extern "C" {
	#include "rebound.h"
}

std::string BasicGravity::updateSim(class Mantle* m, reb_simulation* r) {r->gravity = reb_simulation::REB_GRAVITY_TREE; return "";}
std::string CompensatedGravity::updateSim(class Mantle* m, reb_simulation* r) {r->gravity = reb_simulation::REB_GRAVITY_COMPENSATED; return "";}
std::string TreeGravity::updateSim(class Mantle* m, reb_simulation* r) {
    r->gravity = reb_simulation::REB_GRAVITY_TREE;
    if (opening_angle2.has_value()) {r->opening_angle2 = opening_angle2.value();}
    if (!box.box_size.has_value()) {return "Error: No box size set";}
    reb_simulation_configure_box(r, box.box_size.value(), box.x_boxes.value(), box.y_boxes.value(), box.z_boxes.value());
#ifdef MPI
    reb_mpi_init(r);
#endif
    return "";
}
std::string JacobiGravity::updateSim(class Mantle* m, reb_simulation* r) {r->gravity = reb_simulation::REB_GRAVITY_JACOBI; return "";}
std::string NoGravity::updateSim(class Mantle* m, reb_simulation* r) {r->gravity = reb_simulation::REB_GRAVITY_NONE; return "";}

std::string DirectCollision::updateSim(class Mantle* m, reb_simulation* r) {
    r->collision_resolve = reb_collision_resolve_hardsphere;
    if (collisionContext.has_value()) {m->collisionContext = collisionContext.value();}
    r->collision = reb_simulation::REB_COLLISION_DIRECT;
    return "";
}
std::string LineCollision::updateSim(class Mantle* m, reb_simulation* r) {
    r->collision_resolve = reb_collision_resolve_hardsphere;
    if (collisionContext.has_value()) {m->collisionContext = collisionContext.value();}
    r->collision = reb_simulation::REB_COLLISION_LINE;
    return "";
}
std::string TreeCollision::updateSim(class Mantle* m, reb_simulation* r) {
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->collision = reb_simulation::REB_COLLISION_TREE;
    if (collisionContext.has_value()) {m->collisionContext = collisionContext.value();}
    if (!box.box_size.has_value()) {return "Error: No box size set";}
    reb_simulation_configure_box(r, box.box_size.value(), box.x_boxes.value(), box.y_boxes.value(), box.z_boxes.value());
#ifdef MPI
    reb_mpi_init(r);
#endif
    return "";
}
std::string LineTreeCollision::updateSim(class Mantle* m, reb_simulation* r) {
    r->collision_resolve = reb_collision_resolve_hardsphere;
    r->collision = reb_simulation::REB_COLLISION_LINETREE;
    if (collisionContext.has_value()) {m->collisionContext = collisionContext.value();}
    if (!box.box_size.has_value()) {return "Error: No box size set";}
    reb_simulation_configure_box(r, box.box_size.value(), box.x_boxes.value(), box.y_boxes.value(), box.z_boxes.value());
#ifdef MPI
    reb_mpi_init(r);
#endif
    return "";
}
std::string NoCollision::updateSim(class Mantle* m, reb_simulation* r) {r->collision = reb_simulation::REB_COLLISION_NONE; return "";}

std::string IAS15Integrator::updateSim(class Mantle* m, reb_simulation* r) {
    r->integrator = reb_simulation::REB_INTEGRATOR_IAS15;
    if (epsilon.has_value()) {r->ri_ias15.epsilon = epsilon.value();}
    if (min_dt.has_value()) {r->ri_ias15.min_dt = min_dt.value();}
    if (adaptive_mode.has_value()) {r->ri_ias15.adaptive_mode = adaptive_mode.value();}
    return "";
}
std::string BSIntegrator::updateSim(class Mantle* m, reb_simulation* r) {
    if (!ode_only) {r->integrator = reb_simulation::REB_INTEGRATOR_BS;}
    if (eps_rel.has_value()) {r->ri_bs.eps_rel = eps_rel.value();}
    if (eps_abs.has_value()) {r->ri_bs.eps_abs = eps_abs.value();}
    if (max_dt.has_value()) {r->ri_bs.max_dt = max_dt.value();}
    if (min_dt.has_value()) {r->ri_bs.min_dt = min_dt.value();}
    return "";
}
std::string WHFastIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_WHFAST; return "";}
std::string LeapFrogIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG; return "";}
std::string SEIIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_SEI; return "";}
std::string NoIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_NONE; return "";}

std::string General::updateSim(class Mantle* m, reb_simulation* r) {
    if (G.has_value()) {r->G = G.value();}
    m->MU = MU;
    m->LU = LU;
    m->TU = TU;
    r->testparticle_hidewarnings = silence;
    m->has_ode = has_ode;
    m->has_additional_forces = has_additional_forces;
    m->force_is_velocity_dep = force_is_velocity_dep;
    if (dt.has_value()) {r->dt = dt.value();}
    if (active.has_value()) {r->N_active = active.value();}
    if (softening.has_value()) {r->softening = softening.value();}
    if (removeParticleTest.has_value()) {m->removeParticleTest = removeParticleTest.value();}
    if (simTerminate.has_value()) {m->simTerminate = simTerminate.value();}
#ifdef OPENMP
        if(threads.has_value()) {omp_set_num_threads(threads.value());}
#endif
        return "";
}

