#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <numeric>
#include <any>
#include "mantle.hpp"
#include "mantle_particle.hpp"

enum class WHFastKernel {REB_WHFAST_KERNEL_DEFAULT, REB_WHFAST_KERNEL_MODIFIEDKICK, REB_WHFAST_KERNEL_COMPOSITION, REB_WHFAST_KERNEL_LAZY};
enum class WHFastCoord {REB_WHFAST_COORDINATES_JACOBI, REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC, REB_WHFAST_COORDINATES_WHDS, REB_WHFAST_COORDINATES_BARYCENTRIC};

struct reb_simulation;
struct MConfig { virtual std::string updateSim(class Mantle* m, reb_simulation* r) = 0; virtual ~MConfig() = default;};

struct SimBox {
    std::optional<double> box_size;
    std::optional<int> x_boxes;
    std::optional<int> y_boxes;
    std::optional<int> z_boxes;
    SimBox() {x_boxes = 0.0; y_boxes = 0.0, z_boxes = 0.0;}
};

struct BasicGravity : MConfig { std::string updateSim(class Mantle* m, reb_simulation* r) override; };
struct CompensatedGravity : MConfig { std::string updateSim(class Mantle* m, reb_simulation* r) override; };
struct TreeGravity : MConfig {
    SimBox box;
    std::optional<double> opening_angle2;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct JacobiGravity : MConfig { std::string updateSim(class Mantle* m, reb_simulation* r) override; };
struct NoGravity : MConfig { std::string updateSim(class Mantle* m, reb_simulation* r) override; };

struct DirectCollision : MConfig {
	std::optional<m_coll_context> collisionContext;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct LineCollision : MConfig {
	std::optional<m_coll_context> collisionContext;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct TreeCollision : MConfig {
    SimBox box;
	std::optional<m_coll_context> collisionContext;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct LineTreeCollision : MConfig {
    SimBox box;
	std::optional<m_coll_context> collisionContext;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct NoCollision : MConfig { std::string updateSim(class Mantle* m, reb_simulation* r) override; };

struct IAS15Integrator : MConfig {
    std::optional<double> epsilon;
    std::optional<double> min_dt;
    std::optional<bool> adaptive_mode;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct BSIntegrator : MConfig {
    std::optional<double> eps_rel;
    std::optional<double> eps_abs;
    std::optional<double> min_dt;
    std::optional<double> max_dt;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct WHFastIntegrator : MConfig {
    std::optional<uint32_t> corrector;
    std::optional<bool> corrector2;
    std::optional<WHFastKernel> kernel;
    std::optional<WHFastCoord> coordinates;
    std::optional<bool> recalculate_coordinates_this_timestep;
    std::optional<bool> safe_mode;
    std::optional<bool> keep_unsynchronized;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct LeapFrogIntegrator : MConfig { std::string updateSim(class Mantle* m, reb_simulation* r) override; };
struct SEIIntegrator : MConfig {
    std::optional<double> omega;
    std::optional<double> omegaz;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
};
struct NoIntegrator : MConfig { std::string updateSim(class Mantle* m, reb_simulation* r) override; };

struct General : MConfig {
    std::optional<double> G;
    std::optional<double> MU;
    std::optional<double> LU;
    std::optional<double> TU;
    std::optional<double> dt;
    std::optional<int> active;
    std::optional<double> softening;
    std::optional<int> threads;
    std::optional<bool> has_ode;
    std::optional<bool> has_additional_forces;
    std::string updateSim(class Mantle* m, reb_simulation* r) override;
    General() {has_ode = false; has_additional_forces = false; MU = 1.0; LU = 1.0; TU = 1.0;}
};

#endif // CONFIG_H

