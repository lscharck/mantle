#include "mantle_config.hpp"
extern "C" {
	#include "rebound.h"
}

std::string BasicGravity::updateSim(class Mantle* m, reb_simulation* r) {r->gravity = reb_simulation::REB_GRAVITY_TREE; return "";}
std::string CompensatedGravity::updateSim(class Mantle* m, reb_simulation* r) {r->gravity = reb_simulation::REB_GRAVITY_COMPENSATED; return "";}
std::string TreeGravity::updateSim(class Mantle* m, reb_simulation* r) {
    if (opening_angle2.has_value()) {r->opening_angle2 = opening_angle2.value();}
    if (!box.box_size.has_value()) {
        r->gravity = reb_simulation::REB_GRAVITY_TREE;
        return "Error: No box size set";
    }
    reb_simulation_configure_box(r, box.box_size.value(), box.x_boxes.value(), box.y_boxes.value(), box.z_boxes.value());
    r->gravity = reb_simulation::REB_GRAVITY_TREE;
    return "";
}
std::string JacobiGravity::updateSim(class Mantle* m, reb_simulation* r) {r->gravity = reb_simulation::REB_GRAVITY_JACOBI; return "";}
std::string NoGravity::updateSim(class Mantle* m, reb_simulation* r) {r->gravity = reb_simulation::REB_GRAVITY_NONE; return "";}

std::string DirectCollision::updateSim(class Mantle* m, reb_simulation* r) {
    if (collisionContext.has_value()) {m->collisionContext = collisionContext.value();}
    r->collision = reb_simulation::REB_COLLISION_DIRECT;
    return "";
}
std::string LineCollision::updateSim(class Mantle* m, reb_simulation* r) {
    if (collisionContext.has_value()) {m->collisionContext = collisionContext.value();}
    r->collision = reb_simulation::REB_COLLISION_LINE;
    return "";
}
std::string TreeCollision::updateSim(class Mantle* m, reb_simulation* r) {
    if (collisionContext.has_value()) {m->collisionContext = collisionContext.value();}
    if (!box.box_size.has_value()) {
        r->collision = reb_simulation::REB_COLLISION_TREE;
        return "Error: No box size set";
    }
    reb_simulation_configure_box(r, box.box_size.value(), box.x_boxes.value(), box.y_boxes.value(), box.z_boxes.value());
    r->collision = reb_simulation::REB_COLLISION_TREE;
    return "";
}
std::string LineTreeCollision::updateSim(class Mantle* m, reb_simulation* r) {
    if (collisionContext.has_value()) {m->collisionContext = collisionContext.value();}
    if (!box.box_size.has_value()) {
        r->collision = reb_simulation::REB_COLLISION_LINETREE;
        return "Error: No box size set";
    }
    reb_simulation_configure_box(r, box.box_size.value(), box.x_boxes.value(), box.y_boxes.value(), box.z_boxes.value());
    r->collision = reb_simulation::REB_COLLISION_LINETREE;
    return "";
}
std::string NoCollision::updateSim(class Mantle* m, reb_simulation* r) {r->collision = reb_simulation::REB_COLLISION_NONE; return "";}

std::string IAS15Integrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_IAS15; return "";}
std::string BSIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_BS; return "";}
std::string WHFastIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_WHFAST; return "";}
std::string LeapFrogIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_LEAPFROG; return "";}
std::string SEIIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_SEI; return "";}
std::string NoIntegrator::updateSim(class Mantle* m, reb_simulation* r) {r->integrator = reb_simulation::REB_INTEGRATOR_NONE; return "";}

std::string General::updateSim(class Mantle* m, reb_simulation* r) {
    if(G.has_value()) {r->G = G.value();}
    if(MU.has_value()) {m->MU = MU.value();}
    if(LU.has_value()) {m->LU = LU.value();}
    if(TU.has_value()) {m->TU = TU.value();}
    if(dt.has_value()) {r->dt = dt.value();}
    if(active.has_value()) {r->N_active = active.value();}
    if(softening.has_value()) {r->softening = softening.value();}
    if(has_ode.has_value()) {m->has_ode = has_ode.value();}
    if(has_additional_forces.has_value()) {m->has_additional_forces = has_additional_forces.value();}
#ifdef OPENMP
        if(threads.has_value()) {omp_set_num_threads(threads.value());}
#endif
        return "";
}

